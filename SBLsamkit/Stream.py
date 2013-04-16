####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""Base stream class"""

import os
import re
import sys
import types
from copy import deepcopy
from numpy import inf
from Operators import *
from Utils import *
from Exceptions import *
from Stats import compute_rname_stats

import random

class SAMStream(object):
    """
    Basic SAM Stream class describing common stream/filter operations.  
    Specific streams are designed as subclasses of this basic filter.
    The base SAMStream class allows all recognized LL classes.
    """

    def __init__(self, myname, optionString):
        """
        Sets up default parameters.  No LL, suffix, or operators defined.
        Default 'ss' selector is 'AA'
        """
        self.ops = { 'LL'      : '',	# Location mapping filter
                     'ss'      : SS_AA,	# Strand mapping filter
                     'suffix'  : '',	# Output filename suffix
                     'csvExt'  : '',    # Outut CSV file extension
                     'ops'     : [],	# Operator strings
                     'rnames'  : [],	# List of (RNAME,minPos,maxPos) tuples
                     'RNpairs' : False, # rnames are interpreted as ranges pairs
                     'rngfile' : [],    # Range file or range pair file details
                     'rngops'  : {},    # flank,scope rngf operators
                     'iLrange' : [],    # Single (iLmin,iLmax) tuple
                     'usePCT'  : False, # Flag set if --percentile option used
                     'pctile'  : {},    # Percentile variable dictionary
                     'pctLim'  : {},    # Percentile limits dictionary
                     'sample'  : (0,0), # For sample,N,M operator
                     'filters' : [],    # Special filters
                   }

        self.failsafe = "invalid.sam"
        self.name = myname
        self.optionString = optionString
        self.fileswritten = set()
        self.outputfilenames = set()
        self.count = 0
        self.samplecount = 0
        self.warnings = []
        self.errors = []
        self.stats = { 'sizes'    : [],
                       'rnames'   : {},
                       'RNsingle' : {},
                       'RNpairs'  : {},
                     }
        self.ops['rngops'] = { 'flank'     : 0,
                               'footprint' : RNG_FP_L1,
                               'repeats'   : False,
                             }
        self.globalstats = {}
        for LL in LLnames:
            self.stats[LL] = 0
            self.globalstats[LL] = 0
        for SS in SSnames:
            self.stats[SS.lower()] = 0
            self.globalstats[SS.lower()] = 0
        self.stats['invalid'] = 0


    def op(self, operator):
        """
        Convenience method which returns whether the passed-in operator
        string exists in self.ops['ops'].  If the operator is OP_ORDER ('order')
        it will have a trailing secondary operator (e.g. order:LL).
        """
        if operator is None:
            return True

        if operator == OP_ORDER:
            for op in self.ops['ops']:
                if op.startswith(OP_ORDER) and ':' in op:
                    order,mode = op.split(':')
                    if mode in ORDnames:
                        return True
            return False

        return operator in self.ops['ops']


    def orderOp(self):
        """
        Returns a list of tuples containing each order operator (if any) and
        the sub-stream which it will operator on
        """
        ops = []
        if self.op(OP_ORDER):
            for op in self.ops['ops']:
                if op.startswith(OP_ORDER):
                    ops.append( (self.leftStream(op), op) )
        return ops


    def leftStream(self, operator):
        """
        Return the substream declared to the left of a specific operator.  So
        for instance if an operator sequence looks like "sam,order,LL,fastq"
        and we passed "order:LL" as the operator this should return "sam".  It's
        a way of letting us apply an operator to a single output sub-stream.
        Note that if no sub-stream was declared to the left of the provided
        operator it's assumed to be "sam"
        """
        substreams = (OP_SAM, OP_FASTQ)
        idx = self.ops['ops'].index(operator)
        for i in range(idx,0,-1):
            if self.ops['ops'][i] in substreams:
                return self.ops['ops'][i]
        return OP_SAM
        
        
    def applyEX(self, matchedFlag):
        """
        Checks if the OP_EX operator is in effect, and if so returns the
        inverse of the provided boolean flag, otherwise it returns the same
        value passed in.
        """
        if self.op(OP_EX):
            return not matchedFlag
        else:
            return matchedFlag


    def recordStats(self, pair, matchedFlag):
        """Examinines incoming pair and records counts for LL (all pairs) and 
        SS (UUE pairs only) as well as a list and count of rnames and rname 
        pairs, as appropriate"""

        # Determine LL type
        for LL in LLnames:
            #exec "pairType = pair.%s" % LL
            pairType = self.evaluatePairMapping(pair, LL)
            if pairType:
                self.globalstats[LL] += 1
                if matchedFlag:
                    self.stats[LL] += 1

        # Determine SS type
        for SS in SSnames:
            #exec "pairType = pair.%s" % SS
            pairType = self.evaluatePairMapping(pair, SS)
            if pairType and pair.UUE:
                self.globalstats[SS.lower()] += 1
                if matchedFlag:
                    self.stats[SS.lower()] += 1

        # Determine RNAME and RNAME pair counts (only for matched records).
        # For UUE streams the RNAME pair is always a match so we count it
        # only once
        if matchedFlag:
            if self.ops['LL'] == LL_UUE:
                rnames_as_keys = (pair.first.RNAME, )
            else:
                rnames_as_keys = (pair.first.RNAME, pair.second.RNAME)
            for key in rnames_as_keys:
                if key in self.stats['RNsingle'].keys():
                    self.stats['RNsingle'][key] += 1
                else:
                    self.stats['RNsingle'][key] = 1
            key = pair.RNAMEpair
            if key in self.stats['RNpairs'].keys():
                self.stats['RNpairs'][key] += 1
            else:
                self.stats['RNpairs'][key] = 1

    def evaluatePairMapping(self, pair, attribute):
        """Multi-processor SBLsamkit exposed massive memory accumulation when
        StreamProcessor tasks used an eval() or exec as a way to determine 
        which object property to examine.  This method circumvents the issue
        by exlicitly calling the named attribute.   It will also raise an
        AttributeError if the attribute isn't a valid LL_ or SS_ operator.
        """
        #====================================================================#
        # IMPORTANT.  In reference to the docstring above, YOU MUST MAKE
        # SURE to update this method if you add or remove any LL_ or SS_
        # mapping value in Operators.py
        #====================================================================#

        result = None
        att = attribute
        negate = False
        if attribute.startswith('^'):
            att = attribute[1:]
            negate = True

        # Location mappings (LL)
        if att == LL_UX: result = pair.UX
        elif att == LL_UNM: result = pair.UnM
        elif att == LL_URM: result = pair.UrM
        elif att == LL_URN: result = pair.UrN
        elif att == LL_URS: result = pair.UrS
        elif att == LL_UU: result = pair.UU
        elif att == LL_UUE: result = pair.UUE
        elif att == LL_UUD: result = pair.UUD
        elif att == LL_UA: result = pair.UA
        elif att == LL_XX: result = pair.XX
        elif att == LL_XA: result = pair.XA
        elif att == LL_MX: result = pair.MX
        elif att == LL_MM: result = pair.MM
        elif att == LL_MNM: result = pair.MnM
        elif att == LL_MRM: result = pair.MrM
        elif att == LL_MRN: result = pair.MrN
        elif att == LL_MRS: result = pair.MrS
        elif att == LL_MME: result = pair.MME
        elif att == LL_MMD: result = pair.MMD
        elif att == LL_MA: result = pair.MA
        elif att == LL_NMNM: result = pair.nMnM
        elif att == LL_RMRM: result = pair.rMrM
        elif att == LL_RNRN: result = pair.rNrN
        elif att == LL_RSRS: result = pair.rSrS
        elif att == LL_NMRM: result = pair.nMrM
        elif att == LL_NMRN: result = pair.nMrN
        elif att == LL_NMRS: result = pair.nMrS
        elif att == LL_NMA: result = pair.nMA
        elif att == LL_RMA: result = pair.rMA
        elif att == LL_RNA: result = pair.rNA
        elif att == LL_RSA: result = pair.rSA
        elif att == LL_QC: result = pair.QC
        elif att == LL_CLC1: result = pair.CLC1
        elif att == LL_AA: result = pair.AA

        # Strand mappings (ss)
        elif att == SS_FR: result = pair.fr
        elif att == SS_RF: result = pair.rf
        elif att == SS_SSD: result = pair.ssd
        elif att == SS_FX: result = pair.fx
        elif att == SS_RX: result = pair.rx
        elif att == SS_FF: result = pair.ff
        elif att == SS_RR: result = pair.rr
        elif att == SS_SSE: result = pair.sse
        elif att == SS_AA: result = pair.aa

        # Invalid mapping attribute - raise exception
        if result is None:
            raise AttributeError("%s: Invalid LL or ss test '%s'" % (self.name,
                                                                     attribute))
        # Valid mapping sttribute.  Return match value
        if negate:
            result = not result
        return result


    def match(self, pair):
        """
        Test pair for a given filter.  First it verifies the LL and ss filters
        pass, and if that check succeeds, it will next call the 
        self.matchWithSizes() method.  If that succeeds, it will finally,
        based on the defined operations for the stream (in the self.ops dict), 
        call one of two specific match methods:

        RNAME range list:       self.matchWithRangeList()
        RNAME range pair list:  self.matchWithRangePairList()
        """

        # Test LL and ss operations (or their inverse) succeed:
        matchedLL = self.evaluatePairMapping(pair, self.ops['LL'])
        matchedSS = self.evaluatePairMapping(pair, self.ops['ss'])
        #if '^' in self.ops['LL']:
        #    exec "matchedLL = not pair.%s" % self.ops['LL'].replace('^','')
        #else:
        #    exec "matchedLL = pair.%s" % self.ops['LL']
        #if '^' in self.ops['ss']:
        #    exec "matchedSS = not pair.%s" % self.ops['ss'].replace('^','')
        #else:
        #    exec "matchedSS = pair.%s" % self.ops['ss']
        if not matchedLL or not matchedSS:
            return False

        # Test NMf,RPf etc.... filters if specified
        for filter in self.ops['filters']:
            if filter[0].replace('^','') in (FP_NMF, FP_NMFA):
                if not self.matchWithNMrange(pair, filter):
                    return False
            elif filter[0].replace('^','') in (FP_RPF, FP_RPFA):
                if not self.matchWithRPrange(pair, filter):
                    return False

        # Test iLrange match:
        if len(self.ops['rnames']) == 0:
            # Final filter test:  PHRED scores
            for filter in self.ops['filters']:
                if filter[0].replace('^','') in (FP_PHREDF, FP_PHREDFA):
                    if not self.matchWithPHRED(pair, filter):
                        return False

            return self.matchWithSize(pair)

        # If RNAME range list(s) or range pair list(s) are specified, apply 
        # correct test based on whether the list is part of a range pair list 
        # or not.  Note that we check the iLrange size test first, since if this
        # fails, no RNAME match testing will be performed
        if not self.matchWithSize(pair):
            return False

        # Range or Range pairs
        if self.ops['RNpairs']:
            rnMatch = self.matchWithRangePairList(pair)
        else:
            rnMatch = self.matchWithRangeList(pair)
        if not rnMatch:
            return False

        # Final filter test:  PHRED scores
        for filter in self.ops['filters']:
            if filter[0].replace('^','') in (FP_PHREDF, FP_PHREDFA):
                if not self.matchWithPHRED(pair, filter):
                    return False

        # All tests pass!
        return True
        

    def matchWithRangeList(self, pair):
        """Compares pair against each rname in self.ops['rnames'] and
        returns True if any match is found on either strand.  If a scope
        has been specified (part of the rngf operator) matching is more
        restrictive:

        'L1':      left-most record must match POS if both records match RNAME
        and matching RNAME record must match POS if just one match

        'R2':      right-most record must match POSR2 if both records match 
        RNAME, matching RNAME record must match POSR2 if just one match

        'L1ANDR2': both must match, leftmost to POSL1, right-most to POSR2
        no match if only one record matches the RNAME

        'L1ORR2':  either left must match POSL1, or right must match POSR2.
        For a single RNAME match, either L1 or R2 matches allowed.
        """
        AllowRepeatMapped = self.ops['rngops']['repeats']
        footprint = self.ops['rngops']['footprint']
        flank = self.ops['rngops']['flank']

        r1 = pair.first
        r2 = pair.second

        # Loop through each configured RNAME,minPos,maxPos tuple:
        matchedRNAME = False

        for rname in self.ops['rnames']:
            testedrname = rname

            # Both records match RNAME
            if r1.matchRNAME(rname) and r2.matchRNAME(rname):
                # Skip if either record is repeat-mapped and 'repeat' is not set
                if r1.rM or r2.rM:
                    if not AllowRepeatMapped:
                        continue
                    if r1.rM and not r1.rMM:
                        continue
                    if r2.rM and not r2.rMM:
                        continue
                if footprint == RNG_FP_L1:
                    if pair.left.matchRNAMEandPOS(rname, flank):
                        matchedRNAME = True
                        break
                elif footprint == RNG_FP_R2:
                    if pair.right.matchRNAMEandPOS(rname, flank, RHS=True):
                        matchedRNAME = True
                        break
                elif footprint == RNG_FP_L1ANDR2:
                    if pair.left.matchRNAMEandPOS(rname, flank) and \
                       pair.right.matchRNAMEandPOS(rname, flank, RHS=True):
                        matchedRNAME = True
                        break
                elif footprint == RNG_FP_L1ORR2:
                    if pair.left.matchRNAMEandPOS(rname, flank) or \
                       pair.right.matchRNAMEandPOS(rname, flank, RHS=True):
                        matchedRNAME = True
                        break

            # Only one record matches RNAME
            elif r1.matchRNAME(rname):
                # Skip if matched record is repeat-mapped & 'repeat' is not set
                if r1.matchRNAME(rname):
                    if not r1.rM or (AllowRepeatMapped and r1.rMM):
                        if footprint == RNG_FP_L1:
                            if r1.matchRNAMEandPOS(rname, flank):
                                matchedRNAME = True
                                break
                        elif footprint == RNG_FP_R2:
                            if r1.matchRNAMEandPOS(rname, flank, RHS=True):
                                matchedRNAME = True
                                break
                        elif footprint == RNG_FP_L1ORR2:
                            if r1.matchRNAMEandPOS(rname, flank) or \
                               r1.matchRNAMEandPOS(rname, flank, RHS=True):
                                matchedRNAME = True
                                break

            # Only one record matches RNAME
            elif r2.matchRNAME(rname):
                # Skip if matched record is repeat-mapped & 'repeat' is not set
                if r2.matchRNAME(rname):
                    if not r2.rM or (AllowRepeatMapped and r2.rMM):
                        if footprint == RNG_FP_L1:
                            if r2.matchRNAMEandPOS(rname, flank):
                                matchedRNAME = True
                                break
                        elif footprint == RNG_FP_R2:
                            if r2.matchRNAMEandPOS(rname, flank, RHS=True):
                                matchedRNAME = True
                                break
                        elif footprint == RNG_FP_L1ORR2:
                            if r2.matchRNAMEandPOS(rname, flank) or \
                               r2.matchRNAMEandPOS(rname, flank, RHS=True):
                                matchedRNAME = True
                                break

        if matchedRNAME:
            try:
                self.stats['rnames'][RnameAsKey(testedrname)] += 1
            except KeyError:
                self.stats['rnames'][RnameAsKey(testedrname)] = 1
        return matchedRNAME


    def matchWithRangePairList(self, pair):
        """Compares pair against each pair of rnames in self.ops['rnames'] and
        returns True if a match is found on both strands where each strand
        matches each rname in the pair.
        """
        AllowRepeatMapped = self.ops['rngops']['repeats']
        flank = self.ops['rngops']['flank']

        r1 = pair.first
        r2 = pair.second

        matchedRNAME = False
        for (rname1,rname2) in IterateInPairs(self.ops['rnames']):
            testedrnames = (rname1, rname2)
            if r1.matchRNAMEandPOS(rname1, flank) and r2.matchRNAMEandPOS(rname2, flank):
                if (not r1.rM and not r2.rM) or AllowRepeatMapped:
                    matchedRNAME = True
                    break
            if r1.matchRNAMEandPOS(rname2, flank) and r2.matchRNAMEandPOS(rname1, flank):
                if (not r1.rM and not r2.rM) or AllowRepeatMapped:
                    matchedRNAME = True
                    break
        if matchedRNAME:
            try:
                self.stats['rnames'][RnamePairAsKey(testedrnames)] += 1
            except KeyError:
                self.stats['rnames'][RnamePairAsKey(testedrnames)] = 1
        return matchedRNAME


    def matchWithSize(self, pair):
        """Simplest match method, this returns True if the insert length of
        of the pair falls between the lower and higher size limits (iLrange)
        inclusive the lower limit (i.e. iLmin <= pair.absiL < iLmax). 
        """
        iLmin,iLmax = self.ops['iLrange']
        return iLmin <= pair.absiL < iLmax

    def matchWithPHRED(self, pair, filter):
        """Filter match method, this returns True if the phred score of each
        record (FP_PHREDF) or either record (FP_PHREDFA) falls between the
        specified min and max values provided in the filter (inclusive)
        """
        if filter[0] == FP_PHREDF or filter[0] == '^' + FP_PHREDF:
            match = -2
        elif filter[0] == FP_PHREDFA or filter[0] == '^' + FP_PHREDFA:
            match = -1
        else:
            return False
        if pair.first.matchPHREDvalue(filter[1]): match += 1
        if pair.second.matchPHREDvalue(filter[1]): match += 1
        if '^' in filter[0]:
            return not match >= 0
        else:
            return match >= 0

    def matchWithNMrange(self, pair, filter):
        """Filter match method, this returns True if the NM:i tag is present
        and the value defined falls within self.ops['NMrange'].  Depending
        of whether the 'filters' (self.ops['filters') contains NMf or NMfa, 
        a True value can be returned if either (NMfa) or both (NMf) records 
        match
        NOTE: Beginning with SBLsamkit-4.10 this method has been corrected to
        only apply the test to mapped records.
        """
        if filter[0] != FP_NMF and filter[0] != '^' + FP_NMF and \
           filter[0] != FP_NMFA and filter[0] != '^' + FP_NMFA:
            return False
        if pair.first.mapped:
            p1match = pair.first.matchNMvalue(filter[1])
        else:
            p1match = True
        if pair.second.mapped:
            p2match = pair.second.matchNMvalue(filter[1])
        else:
            p2match = True
        # Match both or just one?
        if filter[0] == FP_NMF or filter[0] == '^' + FP_NMF:
            finalmatch = p1match and p2match
        else:
            finalmatch = p1match or p2match
        if '^' in filter[0]:
            finalmatch = not finalmatch
        return finalmatch

    def matchWithRPrange(self, pair, filter):
        """Filter match method, this returns True if the ZN:i tag is present
        and the value defined falls within self.ops['RPrange'].  Depending
        of whether the 'filters' (self.ops['filters') contains RPf or RPfa, 
        a True value can be returned if either (RPfa) or both (RPf) records 
        match
        """
        if filter[0] == FP_RPF or filter[0] == '^' + FP_RPF:
            match = -2
        elif filter[0] == FP_RPFA or filter[0] == '^' + FP_RPFA:
            match = -1
        else:
            return False
        if pair.first.matchRPvalue(filter[1]): match += 1
        if pair.second.matchRPvalue(filter[1]): match += 1
        if '^' in filter[0]:
            return not match >= 0
        else:
            return match >= 0

    def next(self, pair, options, countOnly=False):
        """
        Returns the pair which came in if the match method returns false,
        others return None, to prevent the pair getting processed by another,
        following filter. If the corresponding converter has the 'pass'
        operator however, it will return the pair which it has output. Other
        converters may also return an unmodified pair, for example the 'strand'
        operator.  The optional countOnly operator can be set to True in which
        case record pairs are analyzed and stats are recorded, but nothing is
        written
        """

        try:
            matchFound = self.match(pair)
            # A match may be artifically rejected if it exceeds a sample counter
            if matchFound and self.op(OP_SAMPLE):
                self.samplecount += 1
                if self.samplecount > self.ops['sample'][0]:
                    matchFound = False
                if self.samplecount == self.ops['sample'][1]:
                    self.samplecount = 0

            matchFound = self.applyEX(matchFound)
            self.recordStats(pair, matchFound)

            # Find and return matches for SP (or countOnly)
            if countOnly or options.nproc == 1:
                if matchFound:
                    self.stats['sizes'].append(pair.absiL)
                    if countOnly:
                        self.count += 1
                    else:
                        for f in self.outputFiles(options):
                            self.write(pair, f)
                    if self.op(OP_PASS):
                        return pair
                    else:
                        return None
                else:
                    return pair
            # Find and return matches for SP (or countOnly)
            else:
                result = {'matched' : False,
                          'name': self.name,
                          'output' : [],
                          'stats' : None,
                         }
                if matchFound:
                    result['output'] = self.outputFiles(options)
                    self.stats['sizes'].append(pair.absiL)
                    result['stats'] = self.stats
                result['matched'] = matchFound
                return result
        except Exception as e:
            self.writefailsafe(pair, options, e)
            self.stats['invalid'] += 1
            return None


    def examinenext(self, pair, options):
        """
        A convenience wrapper for self.next(), this method calls next() with
        the optional countOnly argument set to True.  This processes the stream
        just as it usually does, but without writing any output.
        """
        return self.next(pair, options, countOnly=True)


    def write(self, pair, fileroot):
        """
        Controls writing the SAM and/or FASTQ format output records, after
        applying whatever transformations are necessary.  Saves output
        filename in a 'set' attribute (only one copy of any item added is
        retained), which will be used by the writeInfoFiles() function (in
        the Processor module) after the run has completed.
        """
        self.count += 1

        # FASTQ Output.  Operators which affect output:
        #
        #  OP_NOUT:   Neither .fastq nor .info files written
        #  OP_INFO:   Only .info file is written
        #  OP_NOINFO: Only *.fastq file(s) written
        #  OP_SH:     Writes a single shuffled fastq file
        #  OP_GZ:     Writes *.fastq file(s) gzip-compressed
        if self.op(OP_FASTQ):
            if not self.op(OP_NOUT) and not self.op(OP_INFO):
                fixed, newfileroot = self.ProcessPair(pair, OP_FASTQ,fileroot)
                if self.op(OP_SH):
                    filename1 = '%s.%s' % (newfileroot, 'sh.fastq')
                    filename2 = None
                else:
                    filename1 = '%s.%s' % (newfileroot, '1.fastq')
                    filename2 = '%s.%s' % (newfileroot, '2.fastq')
                fixed.writeFASTQ((filename1, filename2), self.op(OP_GZ))
                self.fileswritten.add('%s.%s' % (newfileroot, 'fastq'))
                self.outputfilenames.add(filename1)
                self.outputfilenames.add(filename2)
        elif self.op(OP_FASTQPP):
            if not self.op(OP_NOUT) and not self.op(OP_INFO):
                fixed, newfileroot = self.ProcessPair(pair, OP_FASTQPP,fileroot)
                if self.op(OP_SH):
                    filename1 = '%s.%s' % (newfileroot, 'sh.fastq')
                    filename2 = None
                else:
                    filename1 = '%s.%s' % (newfileroot, 'pp.1.fastq')
                    filename2 = '%s.%s' % (newfileroot, 'pp.2.fastq')
                fixed.writeFASTQ((filename1, filename2), self.op(OP_GZ))
                self.fileswritten.add('%s.%s' % (newfileroot, 'fastq'))
                self.outputfilenames.add(filename1)
                self.outputfilenames.add(filename2)

        # SAM Output.  Operators which affect output:
        #
        #  OP_NOUT:   Neither .sam nor .info files written
        #  OP_INFO:   Only .info file is written
        #  OP_NOINFO: Only .sam file written
        #  OP_SE:     Writes a single-ended sam file (not supported yet)
        #  OP_GZ:     Writes *.fastq file(s) gzip-compressed
        if self.op(OP_SAM):
            if not self.op(OP_NOUT) and not self.op(OP_INFO):
                fixed, newfileroot = self.ProcessPair(pair, OP_SAM, fileroot)
                filename = '%s.%s' % (newfileroot, 'sam')
                fixed.writeSAM(filename, self.op(OP_GZ))
                self.fileswritten.add(filename)
                self.outputfilenames.add(filename)
        elif self.op(OP_SAMPP):
            if not self.op(OP_NOUT) and not self.op(OP_INFO):
                fixed, newfileroot = self.ProcessPair(pair, OP_SAMPP,fileroot)
                filename = '%s.%s' % (newfileroot, 'pp.sam')
                fixed.writeSAM(filename, self.op(OP_GZ))
                self.fileswritten.add(filename)
                self.outputfilenames.add(filename)


    def writefailsafe(self, pair, options, err):
        """Attempt to write pair to a "error" stream for later analysis without
        terminating the running process."""
        for dir in options.outputDirectories:
            filename = os.path.join(dir, self.failsafe)
            fp = open(filename, "a")
            print >> fp, "# Stream %s: %s" % (self.name, err)
            fp.close()
            pair.writeSAM(filename)

        
    def ProcessPair(self, pair, outputoperator, fileroot):
        """
        Processes output pair to apply proper pair conversions if requested,
        and reverse-complementing strands for FASTQ output etc.
        """
        paircopy = pair.copy()

        # Convert to Proper Pair and/or Reverse Complemented (FASTQ only) 
        if paircopy.UUE:
            if outputoperator == OP_SAMPP:
                # TEMPFIX: trap for stampy records failing PP conversion
                try:
                    paircopy = paircopy.proper_pair(OP_SAM, self.ops)
                except ValueError:
                    fileroot += '.__stampy_PP__'
                except:
                    raise
            elif outputoperator == OP_FASTQPP:
                paircopy = paircopy.proper_pair(OP_FASTQ, self.ops)
            elif outputoperator == OP_FASTQ:
                paircopy = paircopy.reverse_complement()

        # Convert to ordered pair.
        for eachOp in self.orderOp():
            format, op = eachOp
            if self.op(op) and format == outputoperator:
                fileroot = '%s.ord.%s' % (fileroot, op.split(':')[1])
                paircopy = paircopy.ordered_pair(self, op)
                break

        return paircopy, fileroot


    def outputFiles(self, opts):
        """
        Generate a list of one or more output files derived from the
        opts.output_dir and opts.inputfile attributes.  If opts.output_dir
        is None, an output directory named from the current date & time is
        generated
        """
        filename = os.path.split(opts.inputfile)[1]
        fileroot,fileext = os.path.splitext(filename)
        if fileext == '.gz':
            fileroot,fileext = os.path.splitext(fileroot)

        outputfilename = '%s%s' % (fileroot, self.streamSuffix)
        outputfiles = []
        for dir in opts.outputDirectories:
            outputfiles.append(os.path.join(dir, outputfilename))
        return outputfiles


    def parseOperatorTokens(self, tokens):
        """
        The default option string reads follows:
        --stN=suffix,[LL,][ss,][iL,min,max][,rname list|rname pair][,op1,op2,...]
        """
        # (1) Extract suffix, accepts an empty string
        if len(tokens) > 0:
            self.ops['suffix'] = tokens.pop(0)
        if self.ops['suffix'] != '':
            if not self.ops['suffix'].startswith('.'):
                self.ops['suffix'] = '.%s' % self.ops['suffix']

        # (2) Extract LL, ss, and/or filter designator(s) in any order
        while len(tokens) > 0 and (FoundIn(tokens[0], LLnames, True) or 
                                   FoundIn(tokens[0], SSnames, True) or
                                   FoundIn(tokens[0], FPnames, True)):
            if FoundIn(tokens[0], LLnames, True):
                self.ops['LL'] = tokens.pop(0)
            elif FoundIn(tokens[0], SSnames, True):
                self.ops['ss'] = tokens.pop(0)
            elif FoundIn(tokens[0], FPnames, True):
                filter = tokens.pop(0)
                range,tokens = self.extractFilterIntRange(tokens, filter)
                self.ops['filters'].append((filter, range))

        # (3) Extract iLrange
        if len(tokens) > 0:
            if tokens[0] == 'iL':
                tokens.pop(0)
                tokens = self.extractILRange(tokens)

        # (3) Extract RNAMES
        if len(tokens) > 0:
            tokens = self.extractRnames(tokens)

        # SPECIAL CASE:  pp operator is actually a "sub-operator" which may only
        # appear directly *after* an OP_SAM or OP_FASTQ operator.
        corrected_tokens = []
        for tok in tokens:
            if tok == OP_PP:
                try:
                    if corrected_tokens[-1] in (OP_SAM, OP_FASTQ):
                        corrected_tokens[-1] += ',%s' % OP_PP
                    else:
                        raise IndexError
                except IndexError:
                    raise OrderModeParserException("The 'pp' operator may only occur immediately after\na 'sam' or 'fastq' operator", self.name)
            else:
                corrected_tokens.append(tok)
        tokens = corrected_tokens

        # Remaining tokens comprise the operator strings. Any string following
        # the "order" operator should be appended to "order" if itself is not.
        # The "csv" token must be followed by an integer in range 0-99.  The
        # "sample" token must be followed by two integers N and M where N < M
        self.ops['ops'] = []
        index = 0
        i = 0
        while i < len(tokens):
            op = tokens[i]
            if op not in OPnames:
                if len(self.ops['ops']) > 0 and self.ops['ops'][-1] == OP_ORDER:
                    if op in ORDnames:
                        self.ops['ops'][-1] = "%s:%s" % (self.ops['ops'][-1],op)
                        i += 1
                        continue
                    raise OrderModeParserException("Invalid order mode '%s'"%op,
                                                   self.name)
                raise OperatorParserException("Invalid operator '%s'" % op,
                                              self.name)
            # Handle CSV operator
            if op == OP_CSV:
                try:
                    csvnum = tokens[i+1]
                except:
                    raise OperatorParserException(
                                  "Missing N for '%s' operator" % op, self.name)
                try:
                    csvnum = int(csvnum)
                except ValueError:
                    raise OperatorParserException(
                          "Invalid N '%s' for '%s' operator" % (csvnum, op),
                           self.name)
                if not 0 <= csvnum <= 99:
                    raise OperatorParserException(
                          "Invalid N (%d outside range 0-99) for '%s' operator" % (csvnum, op), self.name)
                self.ops['csvExt'] = '%d.%s' % (csvnum, OP_CSV)

                i += 2
                self.ops['ops'].append(op)
                continue

            # Handle SAMPLE operator
            if op == OP_SAMPLE:
                try:
                    sampleN = Integer(tokens[i+1])
                except IndexError:
                    raise OperatorParserException(
                                  "Missing N snd M for '%s' operator" % op, self.name)
                except ValueError:
                    raise OperatorParserException(
                          "Invalid N '%s' for '%s' operator" % (tokens[i+1], op),
                           self.name)
                try:
                    sampleM = Integer(tokens[i+2])
                except IndexError:
                    raise OperatorParserException(
                                  "Missing M for '%s' operator" % op, self.name)
                except ValueError:
                    raise OperatorParserException(
                          "Invalid M '%s' for '%s' operator" % (tokens[i+2], op),
                           self.name)
                if not sampleN < sampleM:
                    raise OperatorParserException(
                          "Invalid N (%d must be < %d) for '%s' operator" % (sampleN, sampleM, op), self.name)
                if sampleN <= 0 or sampleM <= 0:
                    raise OperatorParserException(
                          "Invalid value: M or N is <= 0 for '%s' operator" % op, self.name)
                self.ops['sample'] = (sampleN,sampleM)

                i += 3
                self.ops['ops'].append(op)
                continue

            # Bump token index for next round
            self.ops['ops'].append(op)
            i += 1


        # Verify order mode is valid
        if self.op(OP_ORDER):
            for op in self.ops['ops']:
                if op.startswith(OP_ORDER):
                    order,mode = op.split(':')
                    if mode == ORD_RNG and (len(self.ops['rnames']) == 0 or
                                            self.ops['RNpairs']):
                        raise OrderModeParserException(
                              "Invalid order mode '%s' (no ranges)"% op.replace(':',','), self.name)
                    if mode == ORD_RNGP and (len(self.ops['rnames']) == 0 or
                                             not self.ops['RNpairs']):
                        raise OrderModeParserException(
                              "Invalid order mode '%s' (no range pairs)"% op.replace(':',','), self.name)
                    
        # Verify 'sam' or 'fastq' is specified (NO default)
        if OP_FASTQ not in self.ops['ops'] and \
           OP_FASTQPP not in self.ops['ops'] and \
           OP_SAM not in self.ops['ops'] and \
           OP_SAMPP not in self.ops['ops']:
            raise OrderMissingFormatParserException(
                    "No output format (sam,fast) specified for '%s'"% self.name)

        # Set default iLrange, LL and operators as applicable
        self.setDefaults()


    def extractILRange(self, tokens):
        """
        Extracts iLmin,iLmax values from the operator list.  Supports variable
        substitution (as defined by the --percentile option) and basic
        math operations (+, -, *, /).  Throws an exception if iLmin >= iLmax
        or if the iLrange values and/or expression is not valid.
        """
        iLmin = None
        iLmax = None
        if len(tokens) > 1:
            tok0,tok1 = tokens[:2]
            try:
                iLmin = Integer(tok0)
                iLmax = Integer(tok1)
                if iLmax < 0:
                    iLmax = inf
                if iLmin >= iLmax:
                    raise iLRangeParserException("iLmin >= iLmax: \"%s,%s\"" %
                                                 (tok0, tok1), self.name)
                
                self.ops['iLrange'] = [iLmin, iLmax]
                tokens.pop(0)
                tokens.pop(0)
                return tokens
            except ValueError:
                # If the --percentile option was used, we will continue and
                # see if the iLmin or iLmax string is a valid expression
                if not self.ops['usePCT']:
                    return tokens
                    #raise iLRangeParserException("Illegal iLrange: \"%s,%s\"" %
                    #                             (tok0, tok1), self.name)

            # Avoid regular expression exceptions from RNAME ranges/range pairs
            # getting scanned as an iLrange, if no iLrange was actually 
            # specified
            if '[' in tokens[0] or '[' in tokens[1]:
                return tokens
            try:
                iLmin = MakeExpression(tokens[0])
            except:
                raise iLRangeParserException("Invalid iLrange expression '%s'" %
                                             tokens[0], self.name)
            try:
                iLmax = MakeExpression(tokens[1])
            except:
                raise iLRangeParserException("Invalid iLrange expression '%s'" %
                                             tokens[1], self.name)
            self.ops['iLrange'] = [iLmin, iLmax]
            tokens.pop(0)
            tokens.pop(0)

        return tokens

    def extractFilterIntRange(self, tokens, filter):
        """
        Extracts Nmin,Nmax values from the operator list. Unlike the more
        sophisticated extractILRange() method, this only accepts two integer
        values (represented of course as strings) as a tuple (Nmin,Nmax),
        where Nmax must be >= Nmin.  Throws an exception if the first two
        tokens in the list do not fit these criterea.  Returns the integer range
        as a list, along with the remaining tokens.
        """
        Nmin = None
        Nmax = None
        if len(tokens) > 1:
            tok0,tok1 = tokens[:2]
            try:
                Nmin = Integer(tok0)
                Nmax = Integer(tok1)
            except ValueError:
                # If the --percentile option was used, we will continue and
                # see if the iLmin or iLmax string is a valid expression
                if not self.ops['usePCT']:
                    raise IntRangeParserException(
                                "'%s' specified with an invalid range" % filter,
                                self.name)
                try:
                    Nmin = self.evaluatedFloat(tok0)
                    Nmax = self.evaluatedFloat(tok1)
                except:
                    raise IntRangeParserException(
                                "'%s' specified with an invalid range" % filter,
                                self.name)
            if Nmax < 0:
                Nmax = inf
            if Nmin > Nmax:
                raise IntRangeParserException("Nmin > Nmax: \"%s,%s\"" %
                                              (tok0, tok1), self.name)
                
            tokens.pop(0)
            tokens.pop(0)
            return [Nmin, Nmax], tokens
        else:
            raise IntRangeParserException("NMf or RPf specified without a valid range", self.name)

    def extractRnames(self, tokens):
        """
        Extracts Range lists from the operator list.  A Range is defined as a 
        tuple or one or more (RNAME,minPos,maxPos) tuples, where minPos and 
        maxPos must be integers.  Also accepts range or range pair lists in
        seperate files via the 'rngf' and 'rngpf' operators
        Optional flank,footprint,and repeat subcommands must be extracted 
        also (range, and rangefile only)
        """

        # Characters which are used to denote range lists and range pair lists
        rangeTokens = ( '[', ']',)

        # Check for valid RNAME range argument
        if tokens[0] not in RNG_RANGES and \
           tokens[0] not in RNG_RANGEFILES:
            return tokens

        if tokens[0] == RNG_RNGPAIR or tokens[0] == RNG_RNGPAIRF:
            self.ops['RNpairs'] = True

        # Check for rngf or rngpf (range/range pair file arguments)
        range_token = tokens[0]
        if tokens[0] in RNG_RANGEFILES:
            try:
                rnametokens, rnameHeaders = TokensFromCSV(tokens[1])
                self.ops['rngfile'].append((tokens[1], rnameHeaders))
            except IndexError:
                raise RnameListParserException(
                                   'Missing filename for %s (stream %s)' % (
                                                       tokens[0], self.name))
            except IOError:
                raise RnameListParserException(
                             'Unable to open file "%s" for %s (stream %s)' % (
                                             tokens[1], range_token, self.name))
            except:
                raise RnameListParserException(
                       'Unknown error opening file "%s" for %s (stream %s)' % (
                                              tokens[1], range_token,self.name))
            tokens.pop(0)
            tokens.pop(0)
            # CONVERT FILE CONTENTS INTO TOKEN LIST
            tokens = rnametokens + tokens
        else:
            tokens.pop(0)

        # Beginning at the end of the token list, locate last integer in the 
        # list, searching backwards.  If the last integer is preceded by 
        # the csv or sample operators (or another integer and *then* then
        # sample operator), or the flank subcommand, skip and keep going
        rindex = None
        findex = None
        for i in range(len(tokens)-1, -1, -1):
            try:
                tokenStr = tokens[i]
                for ptoken in rangeTokens:
                    tokenStrRep = tokenStr.replace(ptoken, '')
                Integer(tokenStrRep)
                rindex = i + 1
                if i > 0 and tokens[i-1] == OP_CSV:
                    continue
                if i > 0 and tokens[i-1] == OP_SAMPLE:
                    continue
                if i > 1 and tokens[i-2] == OP_SAMPLE:
                    continue
                if i > 0 and tokens[i-1] == RNG_FLANK:
                    continue
                break
            except ValueError:
                # If the --percentile option was used, we will continue and
                # see if the the tokens are valid expressions
                if not self.ops['usePCT']:
                    continue
                try:
                    self.evaluatedFloat(tokenStr)
                    rindex = i + 1
                    break
                except:
                    continue
        tokenStr1 = tokenStr

        # Now search again from the front to find the first integer.
        for i in range(0,len(tokens)):
            try:
                tokenStr = tokens[i]
                for ptoken in rangeTokens:
                    tokenStrRep = tokenStr.replace(ptoken, '')
                Integer(tokenStrRep)
                findex = i-1
                break
            except ValueError:
                # If the --percentile option was used, we will continue and
                # see if the the tokens are valid expressions
                if not self.ops['usePCT']:
                    continue
                try:
                    self.evaluatedFloat(tokenStr)
                    findex = i - 1
                    break
                except:
                    continue
        tokenStr2 = tokenStr

        # Extract RNAME, minPos, maxPos
        if findex is not None and rindex is not None and findex < rindex:
            # RNAME might itself be an integer so try and compensate for
            # that here
            if (rindex - findex - 1) % 3 == 0:
                findex += 1
            if (rindex - findex) % 3 != 0:
                raise iLRangeParserException("Invalid RNAME range detected - " +
                                             "Declare in " +
                                             "the form RNAME,minPos,maxPos",
                                             self.name)
            for i in range(findex, rindex, 3):
                rnName = tokens.pop(findex)
                minPos = tokens.pop(findex)
                maxPos = tokens.pop(findex)
                for ptoken in rangeTokens:
                    rnName = rnName.replace(ptoken, '')
                    minPos = minPos.replace(ptoken, '')
                    maxPos = maxPos.replace(ptoken, '')
                try:
                    minPosI = Integer(minPos)
                except ValueError:
                    minPosI = self.evaluatedFloat(minPos)
                try:
                    maxPosI = Integer(maxPos)
                except ValueError:
                    maxPosI = self.evaluatedFloat(maxPos)
                rname = [rnName, minPosI, maxPosI]
                self.ops['rnames'].append(rname)
                self.stats['rnames'][RnameAsKey(rname)] = 0

        # Check for range subcommands
        while len(tokens) > 0 and tokens[0] in RNG_COMMANDS:
            if tokens[0] == RNG_REPEAT:
                self.ops['rngops']['repeats'] = True
                tokens.pop(0)
            elif tokens[0] == RNG_FLANK:
                try:
                    self.ops['rngops']['flank'] = Integer(tokens[1])
                    if Integer(tokens[1]) < 0 or \
                       Integer(tokens[1]) > 100000000:
                        raise ValueError
                except ValueError:
                    raise RnameListParserException(
                     'Invalid flank value "%s" for %s (stream %s)' % (
                                          tokens[1],range_token,self.name))
                except IndexError:
                    raise RnameListParserException(
                     'Missing flank value for %s (stream %s)' % (
                                          range_token,self.name))
                tokens.pop(0)
                tokens.pop(0)
            else:
                try:
                    self.ops['rngops']['footprint'] = tokens[1]
                except IndexError:
                    raise RnameListParserException(
                         'Missing footprint value for %s (stream %s)' % (
                                                   range_token,self.name))
                if self.ops['rngops']['footprint'] not in RNG_FOOTPRINTS:
                    raise RnameListParserException(
                       'Invalid footprint value "%s" for %s (stream %s)' % (
                                  self.ops['rngops']['footprint'],
                                  range_token, self.name))
                if range_token in (RNG_RNGPAIR, RNG_RNGPAIRF):
                    self.warnings.append(
                        "WARNING! (%s): '%s' subcommand ignored for range pairs" % (self.name, RNG_FOOTPRINT));
                    print self.warnings[-1]
                tokens.pop(0)
                tokens.pop(0)

        # Return remaining tokens
        return tokens


    def setDefaults(self):

        # If no LL was specified, set one according to the following rules:
        #   iLrange not set, range pair list set:	UU
        #   iLrange not set, range list set:		UA
        #   iLrange is set:				UUE
        #   default (needed?)				AA
        if self.ops['LL'] == '':
            if len(self.ops['rnames']) > 0 and self.ops['RNpairs']:
                self.ops['LL'] = LL_UU
            elif len(self.ops['rnames']) > 0:
                self.ops['LL'] = LL_UA
            elif self.ops['iLrange'] != []:
                self.ops['LL'] = LL_UUE
            else:
                self.ops['LL'] = LL_AA

        # If no iLrange is set, use the 0 to infinity range (all sizes match):
        if self.ops['iLrange'] == []:
            self.ops['iLrange'] = [0, inf]


    def copy(self):
        """
        Returns a deep copy of itself.
        """
        return deepcopy(self)


    @property
    def streamSuffix(self):
        """
        Generates a stream suffix (used to construct the output filename) by
        appending each operator name (where those operators dictate an extra
        suffix).  Valid operators which can be appended to the suffix include:
        OP_EX,
        NOTE:  The OP_SH operator also implies a suffix but in this case it
        replaces the .1.fastq and .2.fastq extensions for FASTQ files with a
        single .sh.fastq (shuffled FASTQ file) - see self.write()
        """
        suffixes = []
        for operator in (OP_EX,):
            if self.op(operator):
                suffixes.append(operator[:3])
        suffix = '.'.join(suffixes)
        if suffix != '':
            suffix = '.' + suffix
        # Special case - a stream whose name is 'n/a' returns the
        # suffix unmodified and without a leading period
        streamSUFFIX = '%s%s' % (self.ops['suffix'], suffix)
        if self.name == 'n/a' and streamSUFFIX[0] == '.':
            streamSUFFIX = streamSUFFIX[1:]
        return streamSUFFIX


    @property
    def evaluatediLrange(self):
        """
        This is a property which returns the iLrange list unaltered if
        both values are ints.  If either value is a string it will contain
        a valid python expression, in which case it will be evaluated
        and the result will be returned.  An easy way to assign these back
        to the iLrange list is:
        
        stream.ops['iLrange'] = stream.evaluatediLrange
        """
        iLmin, iLmax = self.ops['iLrange']
        if type(iLmin) is types.ListType:
            for token in iLmin:
                if token in self.ops['pctile'].keys():
                    index = iLmin.index(token)
                    iLmin.pop(index)
                    iLmin.insert(index, str(self.ops['pctile'][token]))
            try:
                exec "iLminVal = %s" % ''.join(iLmin)
            except (SyntaxError, Exception):
                raise iLRangeParserException("Invalid iLrange expression '%s'" %
                                             ''.join(iLmin), self.name)
            iLmin = iLminVal
        if type(iLmax) is types.ListType:
            for token in iLmax:
                if token in self.ops['pctile'].keys():
                    index = iLmax.index(token)
                    iLmax.pop(index)
                    iLmax.insert(index, str(self.ops['pctile'][token]))
            try:
                exec "iLmaxVal = %s" % ''.join(iLmax)
            except (SyntaxError, Exception):
                raise iLRangeParserException("Invalid iLrange expression '%s'" %
                                             ''.join(iLmax), self.name)
            iLmax = iLmaxVal

        # Some sanity checks
        if iLmin < 0: 
            self.warnings.append(
                    "WARNING! iLmin computed as %d, reset to 0" % iLmin
            )
            print self.warnings[-1]
            iLmin = 0
        if iLmax < 0: 
            if iLmax < -1:
                self.warnings.append(
                    "WARNING! iLmax computed as %d, reset to infinity" % iLmax
                )
                print self.warnings[-1]
            iLmax = inf
        if iLmax <= iLmin:
            self.warnings.append(
               "WARNING! iLmax (%d) is less than iLmin (%d), reset to %d,%d" % (
                   iLmax,iLmin,iLmax,iLmin
               )
            )
            print self.warnings[-1]
            (iLmin, iLmax) = (iLmax, iLmin)

        return [iLmin, iLmax]


    @property
    def OutputFiles(self):
        outputfiles = list(self.outputfilenames)
        outputfiles.sort()
        filenames = []
        for file in outputfiles:
            if file is not None:
                filedir, filename = os.path.split(file)
                if self.op(OP_GZ):
                    filenames.append('%s.gz' % filename)
                else:
                    filenames.append(filename)
        return filenames
        

    def evaluatedFloat(self, floatVal):
        """
        This is a property which returns the floatVal unaltered if it is
        an float (or typcasts up if an int).  If the value is a list it 
        should contain a valid python expression, in which case it will be 
        evaluated and the result will be returned. A string is recursively
        passed as a single-item list
        """
        if type(floatVal) is types.FloatType:
            return floatVal
        elif type(floatVal) in (types.IntType, types.LongType):
            return float(floatVal)
        elif type(floatVal) is types.StringType:
            return self.evaluatedFloat([floatVal])
        elif type(floatVal) is types.ListType:
            for token in floatVal:
                if token in self.ops['pctile'].keys():
                    index = floatVal.index(token)
                    floatVal.pop(index)
                    floatVal.insert(index, str(self.ops['pctile'][token]))
            try:
                exec "fVal = %s" % ''.join(floatVal)
            except (SyntaxError, Exception):
                raise ValueError("Invalid float expression '%s'" %
                                             ''.join(floatVal), self.name)
            return fVal
        else:
            raise ValueError("Invalid float expression '%s'" %
                                         ''.join(floatVal), self.name)
    

def ConfigureStreams(opts):
    """
    Applies some configurations to each configured stream (e.g. setting the
    stream.ops['pctile'] dict correctly
    """
    try:
        opts = ProcessPercentileOption(opts)
    except:
        ErrorHandler()
    uniqueNames = []
    for stream in opts.orderedStreams:
        try:
            if stream.name in uniqueNames:
                raise DuplicateStreamException("Stream '%s' " % stream.name +
                                               "declared multiple times")
        except:
            ErrorHandler()
        uniqueNames.append(stream.name)
        if opts.pct:
            stream.ops['pctile'] = opts.pct.copy()
            stream.ops['pctLim'] = opts.pctLim.copy()
            stream.ops['usePCT'] = True
        tokens = Tokenize(stream.optionString)
        try:
            stream.parseOperatorTokens(tokens)
        except:
            ErrorHandler()
