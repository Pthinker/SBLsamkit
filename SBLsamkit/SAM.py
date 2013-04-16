####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""Classes and functions related to SAM Pairs"""

import sys
import re
import os.path
import gzip
import types
from SBLhtseq import *
from Operators import ORD_LL, ORD_SS, ORD_POS, ORD_RNG, ORD_RNGP, OP_SAM, OP_FASTQ
from Utils import IterateInPairs

def ReadSAMFile(filename):
    """
    This function wraps the internals of reading the SAM records, using
    both original HTSeq and modified SBLsamkit.SBLhtseq classes and 
    functions, presenting a simple interface to the caller.  It acts as a 
    simple iterator, yielding a single SAMPair on each next() call.  So for 
    example the caller can simply do:
 
    from SBLsamkit.SAM import ReadSAMFile
    for pair in ReadSAMFile(filename): ...

    """
    sam_reader = SAM_Reader(filename)
    for read1, read2 in pair_SAM_alignments(sam_reader):
        yield SAMPair(read1, read2)


def makeSAMpairFromStringTuple(readpair, reorder=True):
    """
    Convert two strings (assumed to be read in from a SAM file) into a valid
    SAM pair.  Not only does each string get converted into a SAM alignment
    object, but if the optional reorder parameter is True (as it is by default)
    the order of two objects must be adjusted so that whichever has a pe_which 
    attribute of 'first' has to be first in the pair.  See the 
    pair_SAM_alignments() function in HTSeq's __init__.py for more information,
    specificially the sub-function called "proces_paired_reads".  Of course
    this implies only PE reads are aceptable for SBLsamkit as it currently 
    stands.
    """
    rec1 = SAM_Alignment(readpair[0])
    rec2 = SAM_Alignment(readpair[1])
    if reorder and rec1.pe_which == 'second':
        tmp = rec1
        rec1 = rec2
        rec2 = tmp
    return SAMPair(rec1, rec2)


class SAMPair(object):
    """
    A class describing a read pair.  Both reads of a pair are stored, and
    some pair-specific attributes are set, such as the computed distance 
    between their start positions.
    """
    def __init__(self, first, second):

        self.first    = first
        self.second   = second
        self.read     = [self.first, self.second]
        self.strand   = [self.first.STRAND, self.second.STRAND]
        self.position = [first.POS, second.POS]
        self.absiL    = abs(first.POS - second.POS)

        if self.strand == ['+','-'] or self.strand == ['-','+']:
            if self.strand[0] == '+':
                self.posF = first.POS
                self.posR = second.POS
            else:
                self.posR = first.POS
                self.posF = second.POS
        else:
            # FIXME: Verify this is correct for unaligned???
            self.posF = first.POS
            self.posR = second.POS

        # FIXME: This fixes negative sizes in RJP... is it correct??
        # self.iL = abs(self.posR - self.posF)
        # Regression - for now reverting to 0.47.alpha3 behaviour
        self.iL = self.posR - self.posF


    def __str__(self):
        """
        String representation of the SAM pair 
        """
        return "\n" + repr(self.first) + "\n" + repr(self.second)


    def writeSAM(self, outputfile, gzipFlag=False, closeWhenDone=True):
        """
        Basic output method - overridden in the Operators subclasses
        to add proper-pair handling.  
        """
        fp = self.__open__(outputfile, 'a', gzipFlag)
        
        print >> fp, self.first.STRING
        print >> fp, self.second.STRING
        if closeWhenDone:
            fp.close()


    def reverse_complement(self):
        """
        Returns a copy of the self with the appropriate strands 
        reverse-complemented for output.  Note that we assign the
        altered reads to the _read attribute, which is then always assumed
        to contain the sequence in "ready for output" orientation.
        """
        pair = self.copy()
        rec1,rec2 = pair.read
        if rec1.STRAND == '-':
            if not self.first.nM and not self.first.rM:
                pair.first._read = rec1._read.get_reverse_complement()
        if rec2.STRAND == '-':
            if not self.second.nM and not self.second.rM:
                pair.second._read = rec2._read.get_reverse_complement()
        return pair


    def proper_pair(self, target, opts):
        """
        Returns a copy of the self with the appropriate strands modified
        into a virtual proper pair.  If self is not UUE, return an unmodified
        copy
        """
        pair = self.copy()
        if not pair.UUE:
            return pair

        if target == OP_SAM:
            # Check if already proper pair
            if self.first.__getflagbit__(SAM_PROPERPAIR) and \
               self.second.__getflagbit__(SAM_PROPERPAIR):
                return pair
            # Check if STP or CJP:
            if self.fr or self.rf:
                if (self.first.__getflagbit__(SAM_QUERYSTRAND) and not 
                    self.first.__getflagbit__(SAM_MATESTRAND)) or \
                   (self.first.__getflagbit__(SAM_MATESTRAND) and not
                   self.first.__getflagbit__(SAM_QUERYSTRAND)):
                    # Swap if CJP:
                    if self.posF > self.posR:
                        pair.first = self.second.copy()
                        pair.second = self.first.copy()
                        pair.first.__setflagbit__(SAM_QUERYSTRAND, 'toggle')
                        pair.first.__setflagbit__(SAM_MATESTRAND, 'toggle')
                        pair.second.__setflagbit__(SAM_QUERYSTRAND, 'toggle')
                        pair.second.__setflagbit__(SAM_MATESTRAND, 'toggle')
                    uq = 0
                    for r in pair.read:
                        for t in r.TAGS:
                            if t[:2] == ['UQ', 'i']:
                                uq += int(t[2])
                    pptags = [ ['PQ', 'i', '%d' % uq],
                               ['SM', 'i', '151'],
                               ['AM', 'i', '0'],
                             ]
                    pair.first._tags.extend(pptags)
                    pair.second._tags.extend(pptags)
                    pair.first.__setflagbit__(SAM_PROPERPAIR, 'on')
                    pair.second.__setflagbit__(SAM_PROPERPAIR, 'on')
                    return pair
            # Check for RJP
            if self.ff or self.rr:
                if (self.first.__getflagbit__(SAM_QUERYSTRAND) and 
                    self.first.__getflagbit__(SAM_MATESTRAND)) or \
                   (not self.first.__getflagbit__(SAM_MATESTRAND) and not
                   self.first.__getflagbit__(SAM_QUERYSTRAND)):
                    if self.first.POS <= self.second.POS:
                        pair.first.__setflagbit__(SAM_QUERYSTRAND, 'off')
                        pair.first.__setflagbit__(SAM_MATESTRAND, 'on')
                        pair.second.__setflagbit__(SAM_QUERYSTRAND, 'on')
                        pair.second.__setflagbit__(SAM_MATESTRAND, 'off')
                    else:
                        pair.first.__setflagbit__(SAM_QUERYSTRAND, 'on')
                        pair.first.__setflagbit__(SAM_MATESTRAND, 'off')
                        pair.second.__setflagbit__(SAM_QUERYSTRAND, 'off')
                        pair.second.__setflagbit__(SAM_MATESTRAND, 'on')
                    uq = 0
                    for r in pair.read:
                        for t in r.TAGS:
                            if t[:2] == ['UQ', 'i']:
                                uq += int(t[2])
                    pptags = [ ['PQ', 'i', '%d' % uq],
                               ['SM', 'i', '151'],
                               ['AM', 'i', '0'],
                             ]
                    pair.first._tags.extend(pptags)
                    pair.second._tags.extend(pptags)
                    pair.first.__setflagbit__(SAM_PROPERPAIR, 'on')
                    pair.second.__setflagbit__(SAM_PROPERPAIR, 'on')
                    return pair
            raise ValueError("Error Converting to PP (SAM)")
        elif target == OP_FASTQ:
            rc = pair.reverse_complement()
            r1 = rc.first
            r2 = rc.second
            #if self.fr: print "F-R: (%s) %9d,%9d" % (rc.strand,r1.POS,r2.POS),
            #if self.rf: print "R-F: (%s) %9d,%9d" % (rc.strand,r1.POS,r2.POS),
            #if self.ff: print "F-F: (%s) %9d,%9d" % (rc.strand,r1.POS,r2.POS),
            #if self.rr: print "R-R: (%s) %9d,%9d" % (rc.strand,r1.POS,r2.POS),
            #print "rec1 > rec2  %5s  " % (r1.POS > r2.POS), 
            #print "posF > posR  %5s  " % (rc.posF > rc.posR),
            iLmin,iLmax = opts['iLrange']
            # "CJP"
            if self.fr or self.rf:
                if self.posF > self.posR:
                    #print "CJP",
                    rc.first._read =  r1._read.get_reverse_complement()
                    rc.second._read = r2._read.get_reverse_complement()
            # "RJP"
            if self.ff or self.rr:
                #print "RJP",
                if self.rr:
                    #print "(rr)",
                    r1._read = r1._read.get_reverse_complement()
                    r2._read = r2._read.get_reverse_complement()
                if r1.POS < r2.POS:
                    #print "(1<2)"
                    rc.first._read  = r1._read
                    rc.second._read = r2._read.get_reverse_complement()
                elif r1.POS > r2.POS:
                    #print "(2>1)"
                    rc.first._read = r1._read.get_reverse_complement()
                    rc.second._read  = r2._read
            #print
            return rc
     

    def proper_pairOLD(self, target, opts):
        """
        Returns a copy of the self with the appropriate strands modified
        into a virtual proper pair.  If self is not UUE, return an unmodified
        copy
        """
        if target == 'sam':
            if self.first.__getflagbit__(SAM_PROPERPAIR) and \
               self.second.__getflagbit__(SAM_PROPERPAIR):
                return pair

            if self.ff or self.rr:
                if self.posF <= self.posR:
                    pair.first.strand  = '+'
                    pair.second.strand = '-'
                    pair.first.__setflagbit__(SAM_QUERYSTRAND, 'off')
                    pair.first.__setflagbit__(SAM_MATESTRAND, 'on')
                    pair.second.__setflagbit__(SAM_QUERYSTRAND, 'on')
                    pair.second.__setflagbit__(SAM_MATESTRAND, 'off')
                else:
                    pair.first.strand  = '-'
                    pair.second.strand = '+'
                    pair.first.__setflagbit__(SAM_QUERYSTRAND, 'on')
                    pair.first.__setflagbit__(SAM_MATESTRAND, 'off')
                    pair.second.__setflagbit__(SAM_QUERYSTRAND, 'off')
                    pair.second.__setflagbit__(SAM_MATESTRAND, 'on')
            else:
                if self.posF > self.posR:
                    pair.first.__setflagbit__(SAM_QUERYSTRAND, 'toggle')
                    pair.first.__setflagbit__(SAM_MATESTRAND, 'toggle')
                    pair.second.__setflagbit__(SAM_QUERYSTRAND, 'toggle')
                    pair.second.__setflagbit__(SAM_MATESTRAND, 'toggle')
                pair.first.__setflagbit__(SAM_MAPPED, 'on')
                pair.second.__setflagbit__(SAM_MAPPED, 'on')
                uq = 0
                for r in pair.read:
                    for t in r.TAGS:
                        if t[:2] == ['UQ', 'i']:
                            uq += int(t[2])
                pptags = [ ['PQ', 'i', '%d' % uq],
                           ['SM', 'i', '151'],
                           ['AM', 'i', '0'],
                         ]
                pair.first._tags.extend(pptags)
                pair.second._tags.extend(pptags)

        elif target == 'fastq':
            r1,r2 = self.read
            iLmin,iLmax = opts['iLrange']
            # "CJP"
            if self.fr or self.rf:
                if self.posF > self.posR:
                    if iLmin <= -self.iL <= iLmax:
                        if self.fr:
                            pair.first._read = r1._read.get_reverse_complement()
                        if self.rf:
                            pair.second._read =r2._read.get_reverse_complement()
            # "RJP"
            if self.ff or self.rr:
                if iLmin <= abs(self.iL) <= iLmax:
                    if r1.POS < r2.POS:
                        pair.second._read = r2._read.get_reverse_complement()
                    elif r1.POS > r2.POS:
                        pair.first._read = r1._read.get_reverse_complement()

        return pair


    def ordered_pair(self, stream, op):
        """
        Returns a copy of the self with the appropriate strands ordered 
        depending on the order operator mode and stream type. 
        """
        order,mode = op.split(':')
        newpair = self.copy()

        # order,LL:
        # modifies the order based on the LL value.  If the LL type is 
        # asymmetric (not the same in both records of a record pair)  (whether 
        # or not there is an LL field in the stream string). Rules for LL 
        # ordering:
        #  - if only one L value of the pair is U  (e.g. UnM, UrM ), then the U 
        #    record is the A record, 
        #  - else if only one of the records is nM (e.g. nMrM) then the nM 
        #    record is the A record,
	#  - else leave record pair unchanged.
        if mode == ORD_LL:
            if self.UnM or self.UrM or self.UX:
                if self.second.U:
                    newpair = self.swap()
            elif self.nMrM:
                if self.second.nM:
                    newpair = self.swap()
        # order,ss:
        # If the ss type  is asymmetrical, whether or not an ss field is 
        # present in the stream definition string, process  the pair as follows:
        #  - If only one s value of the pair is f  (e.g. fr, fx) then process 
        #    the record such that the f record is the A record
	#  - Else if only one s value of a pair is r (e.g. rx) then process 
        #    the record pair such that the r record is the A record,
	#  - Else leave the record pair unchanged.
        if mode == ORD_SS:
            if self.fr or self.fx:
                if self.second.f:
                    newpair = self.swap()
            elif self.rf or self.rx:
                if self.second.r:
                    newpair = self.swap()
        # order,pos:
        # If the record pair type is UUE and if the pos values of the two 
        # records are different, process the record pair such that the record 
        # with the lower pos, is the A record. Otherwise no modification is made
        if mode == ORD_POS:
            if self.UUE and self.first.POS != self.second.POS:
                if self.first.POS > self.second.POS:
                    newpair = self.swap()
        # order,rng:
        # This string is only valid if a range list with one or more ranges is 
        # part of the stream script.  In such a case, if at least one of the 
        # records is mapped, and is mapped within one of the ranges in the list,
        # and only one of the records meets this criterion, then the record pair
        # is processed such that the range matching record is the A record.  
        # Otherwise no modification is made

        if mode == ORD_RNG:
            # Get flank for determining position match
            flank = stream.ops['rngops']['flank']
            r1 = self.first
            r2 = self.second
            matchedRNAME1 = False
            matchedRNAME2 = False
            for rname in stream.ops['rnames']:
                if r1.matchRNAMEandPOS(rname, flank):
                    matchedRNAME1 = True
                if r2.matchRNAMEandPOS(rname, flank):
                    matchedRNAME2 = True
            if matchedRNAME1 != matchedRNAME2:
                if matchedRNAME2:
                    newpair = self.swap()
        # order,rngp:
        # This string is only valid if a range pair list with one or more pairs
        # is part of the stream script.  In such a case, if the two records in 
        # the pair are mapped to a range pair, then the record pair is processed
        # such that A record is the record which matches the first range in a 
        # range pair definition , and the record which matches the second range
        # in a pair is the B record. Otherwise no modification is made.
        if mode == ORD_RNGP:
            # Get flank for determining position match
            flank = stream.ops['rngops']['flank']
            r1 = self.first
            r2 = self.second
            matchedRNAME = False
            for (rname1,rname2) in IterateInPairs(stream.ops['rnames']):
                testedrnames = (rname1, rname2)
                if r1.matchRNAMEandPOS(rname2, flank) and r2.matchRNAMEandPOS(rname1, flank):
                    newpair = self.swap()
        return newpair


    def writeFASTQ(self, outputfiles, gzipFlag=False, closeWhenDone=True):
        """
        Basic output method - overridden in the Operators subclasses
        to add proper-pair handling.  NOTE: This basic method does NO
        reverse-complement processing.  If outputfiles[1] is None, both
        records will be written to the same file (outputfiles[0]), i.e. a
        shuffled output file.
        """
        fp1 = self.__open__(outputfiles[0], 'a', gzipFlag)
        if outputfiles[1] is None:
            fp2 = fp1
        else:
            fp2 = self.__open__(outputfiles[1], 'a', gzipFlag)

        self.__write_fastq__(self.read[0], fp1, "_1")
        self.__write_fastq__(self.read[1], fp2, "_2")
        if closeWhenDone:
            fp1.close()
            fp2.close()
 

    def __open__(self, filename, flags, gzipped=False):
        """
        Opens a file handle with the specified flags.  If gzipped is True,
        the filename will have '.gz' appended to it and gzip.open() will be
        used instead of open().  In either case, returns the open file handle
        For maximum portabililty, the filename may also be an open filehandle
        (or filehandle-like object such as a StringIO object) in which case 
        it is returned unchanged (and in that case the flags and optional 
        gzipped argument are ignored)
        """
        if type(filename) in (types.FileType, types.InstanceType):
            if hasattr(filename, 'write'):
                return filename

        if gzipped:
            return gzip.open('%s.gz' % filename, flags)
        else:
            return open(filename, flags)


    def __write_fastq__(self, record, filehandle, append_to_name=""):
        """
        This private method handles writing a single SAM record to an
        open file handle in FASTQ format
        """
        fastq_name = record.QNAME
        print >> filehandle, "@" + fastq_name + append_to_name
        print >> filehandle, record.SEQ
        print >> filehandle, "+"
        print >> filehandle, record.QUAL


    def copy(self):
        return SAMPair(self.first.copy(), self.second.copy())


    def swap(self):
        """Returns a copy if self where the strands are swapped"""
        return SAMPair(self.second.copy(), self.first.copy())


    @property
    def RNAMEpair(self):
        """Returns a string containing the sorted RNAME pair comma-separated,
        useful as a dictionary key for example.  The pair of RNAMEs is sorted
        such that the pair will be the same regardless of the order they appear
        in the record pair
        """
        rnames = [self.first.RNAME, self.second.RNAME]
        rnames.sort()
        return ','.join(rnames)


    #-- LL properties ---------------------------------------------------------#
    # Properties for easy query of record location type (LL) of record pairs. 
    # These are generally defined in the same order as they are described in 
    # the Samkit 3 spec document, apart from some internal properties which 
    # are defined last
    @property
    def UX(self):
        """
        This property returns True if one record is U and the other is X
        """
        if self.first.U and self.second.X: return True
        if self.first.X and self.second.U: return True
        return False

    @property
    def UnM(self):
        """
        This property returns True when one record is U and the other is nM
        (in either order)
        """
        if self.first.U  and self.second.nM:  return True
        if self.first.nM and self.second.U:   return True
        return False

    @property
    def UrM(self):
        """
        This property returns True when one record is U and the other is rM
        (in either order)
        """
        if self.first.U  and self.second.rM:  return True
        if self.first.rM and self.second.U:   return True
        return False

    @property
    def UrN(self):
        """
        This property returns True when one record is U and the other is rN
        (in either order); (rN means mapped in novoalign style)
        """
        if self.first.U  and self.second.rN:  return True
        if self.first.rN and self.second.U:   return True
        return False

    @property
    def UrS(self):
        """
        This property returns True when one record is U and the other is rS
        (in either order); (rS means mapped in novoalign style)
        """
        if self.first.U  and self.second.rS:  return True
        if self.first.rS and self.second.U:   return True
        return False

    @property
    def UU(self):
        """
        This property returns True only if both records in the pair are 
        uniquely mapped
        """
        return self.first.U and self.second.U

    @property
    def UUE(self):
        """
        A more specific version of the UU property, where both records must
        be mapped to the same RNAMEs
        """
        if self.UU:
            if self.first.RNAME == self.second.RNAME: 
                return True
            if self.first.RNAME != '=' and self.second.RNAME == '=':
                return True
        return False

    @property
    def UUD(self):
        """
        In a similar vein to UUE, both records must be uniquely mapped, but
        to *different* RNAMEs
        """
        if self.UU:
            if self.first.RNAME != '=' and self.second.RNAME == '=':
                return False
            if self.first.RNAME != self.second.RNAME:
                return True
        return False

    @property
    def UA(self):
        """
        At least one record in the pair is U
        """
        return self.first.U or self.second.U

    @property
    def MX(self):
        """
        One record is M, the other is X, in either order
        """
        if self.first.M and self.second.X: return True
        if self.first.X and self.second.M: return True
        return False

    @property
    def MnM(self):
        """
        One record is M, the other is nM, in either order
        """
        if self.first.M and self.second.nM: return True
        if self.first.nM and self.second.M: return True
        return False

    @property
    def MrM(self):
        """
        One record is M, the other is rM, in either order
        """
        if self.first.M and self.second.rM: return True
        if self.first.rM and self.second.M: return True
        return False

    @property
    def MrN(self):
        """
        One record is M, the other is rN (repeat mapped/novoalign), in either order
        """
        if self.first.M and self.second.rN: return True
        if self.first.rN and self.second.M: return True
        return False

    @property
    def MrS(self):
        """
        One record is M, the other is rS (repeat mapped, stampy), in either order
        """
        if self.first.M and self.second.rS: return True
        if self.first.rS and self.second.M: return True
        return False

    @property
    def MM(self):
        """
        Both records are M
        """
        return self.first.M and self.second.M

    @property
    def MME(self):
        """
        A more specific version of the MM property, where both records must
        be mapped to the same RNAMEs
        """
        if self.MM:
            if self.first.RNAME == self.second.RNAME: 
                return True
            if self.first.RNAME != '=' and self.second.RNAME == '=':
                return True
        return False

    @property
    def MMD(self):
        """
        In a similar vein to MME, both records must be mapped, but
        to *different* RNAMEs
        """
        if self.MM:
            if self.first.RNAME != '=' and self.second.RNAME == '=':
                return False
            if self.first.RNAME != self.second.RNAME:
                return True
        return False

    @property
    def MA(self):
        """
        At least one record in the pair is M
        """
        return self.first.M or self.second.M

    @property
    def nMnM(self):
        """
        Convenience property - matches any record pair where both records are
        not mapped (nM)
        """
        return self.first.nM and self.second.nM

    @property
    def rMrM(self):
        """
        Convenience property - matches any record pair where both records are
        repeat mapped (rM)
        """
        return self.first.rM and self.second.rM

    @property
    def rNrN(self):
        """
        Convenience property - matches any record pair where both records are
        repeat mapped in novoalign parlance (rN)
        """
        return self.first.rN and self.second.rN

    @property
    def rSrS(self):
        """
        Convenience property - matches any record pair where both records are
        repeat mapped in stampy parlance (rS)
        """
        return self.first.rS and self.second.rS

    @property
    def nMrM(self):
        """
        Convenience property - matches any record pair where one record is
        unmapped (nM) and the other is repeat mapped (rM), in either order
        """
        if self.first.rM and self.second.nM: return True
        if self.first.nM and self.second.rM: return True
        return False

    @property
    def nMrN(self):
        """
        Convenience property - matches any record pair where one record is
        unmapped (nM) and the other is repeat mapped in novoalign paralace (rN), 
        in either order
        """
        if self.first.rN and self.second.nM: return True
        if self.first.nM and self.second.rN: return True
        return False

    @property
    def nMrS(self):
        """
        Convenience property - matches any record pair where one record is
        unmapped (nM) and the other is repeat mapped in stampy parlance (rS), 
        in either order
        """
        if self.first.rS and self.second.nM: return True
        if self.first.nM and self.second.rS: return True
        return False

    @property
    def nMA(self):
        """
        Convenience property - matches any record pair where at least one
        record is unmapped (nM)
        """
        return self.first.nM or self.second.nM

    @property
    def rMA(self):
        """
        Convenience property - matches any record pair where at least one
        record is repeat mapped (rM)
        """
        return self.rNA or self.rSA

    @property
    def rNA(self):
        """
        Convenience property - matches any record pair where at least one
        record is repeat mapped in novoalign parlance (rN)
        """
        return self.first.rN or self.second.rN

    @property
    def rSA(self):
        """
        Convenience property - matches any record pair where at least one
        record is repeat mapped in stampy parlance (rS)
        """
        return self.first.rS or self.second.rS

    @property
    def XX(self):
        """
        Convenience property - matches any record pair where *neither* record
        is U
        """
        return self.first.X and self.second.X

    @property
    def XA(self):
        """
        Convenience property - matches any record pair where at least one 
        record is X (i.e. either nM or rM)
        """
        return self.first.X or self.second.X

    @property
    def AA(self):
        """
        Convenience property - matches any record pair, since each record will
        always return True for the A property
        """
        if self.first.A or self.second.A:
            return True
        return False
        #return self.first.A or self.second.A


    @property
    def QC(self):
        """
        Returns true if either record in the pair failed any QC tests
        """
        return self.first.QC or self.second.QC


    @property
    def CLC1(self):
        """
        One or both records contain an 'I' in the CIGAR string
        """
        return self.first.CLC1 or self.second.CLC1

    @property
    def UQC(self):
        """
        ** Internally used ONLY **
        This property returns True when one record is U and the other is tagged
        as a "QC" failed read (in either order)
        """
        if self.first.U  and self.second.QC:  return True
        if self.first.QC and self.second.U:   return True
        return False

    @property
    def UCLC1(self):
        """
        ** Internally used ONLY **
        This property returns True when one record is U and the other is tagged
        as a "CLC1" read (in either order)
        """
        if self.first.U    and self.second.CLC1:  return True
        if self.first.CLC1 and self.second.U:   return True
        return False


    #-- SS properties ---------------------------------------------------------#
    # Properties for easy query of strand mapping (ss) of record pairs. 
    # These are generally defined in the same order as they are described in 
    # the Samkit 4 spec document, apart from some internal properties which 
    # are defined last.  
    @property
    def fr(self):
        """
        This property returns True ONLY if the leftmost strand in the pair is
        forward AND the rightmost strand is reverse (+/-)
        """
        if self.first.POS <= self.second.POS:
        	return self.MME and (self.first.f and self.second.r)
	else:
        	return self.MME and (self.first.r and self.second.f)

    @property
    def rf(self):
        """
        This property returns True ONLY if the leftmost strand in the pair is
        reverse AND the rightmost strand is forward: (-/+)
        """
        if self.second.POS < self.first.POS:
        	return self.MME and (self.first.f and self.second.r)
	else:
        	return self.MME and (self.first.r and self.second.f)

    @property
    def ssd(self):
        """
        This property returns True if the pair is either fr or rf 
        """
        return self.MME and (self.fr or self.rf)

    @property
    def fx(self):
        """
        This property returns True if the U strand in a UX pair is mapped
        to the forward strand
        """
        if self.MX:
            if self.first.M:
                return self.first.f
            else:
                return self.second.f
        return False

    @property
    def rx(self):
        """
        This property returns True if the U strand in a UX pair is mapped
        to the forward strand
        """
        if self.MX:
            if self.first.M:
                return self.first.r
            else:
                return self.second.r
        return False

    @property
    def ff(self):
        """
        This property indicates whether both strands are forward-facing (+/+)
        """
        return self.MME and (self.first.f and self.second.f)

    @property
    def rr(self):
        """
        This property indicates whether both strands are reverse-facing (-/-)
        """
        return self.MME and (self.first.r and self.second.r)

    @property
    def sse(self):
        """
        This property returns True if the pair is either ff or rr 
        """
        return self.MME and (self.ff or self.rr)

    @property
    def aa(self):
        """
        This property returns True always 
        """
        return True

    @property
    def left(self):
        """Returns the left-most record of the pair"""
        if self.first.POS <= self.second.POS:
            return self.first
        else:
            return self.second

    @property
    def right(self):
        """Returns the right-most record of the pair"""
        if self.first.POS <= self.second.POS:
            return self.second
        else:
            return self.first
