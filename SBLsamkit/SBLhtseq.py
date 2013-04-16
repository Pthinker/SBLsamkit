####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
Custom SAM_Alignment classes and functions which enhance the functionality
of the HTSeq python package
"""

import re
import types
import HTSeq
import warnings
from numpy import inf

SAM_NONE          = 0x0000	# null
SAM_MAPPED        = 0x0001	# read is paired in sequencing
SAM_PROPERPAIR    = 0x0002	# read is mapped in proper pair
SAM_UNMAPPED      = 0x0004	# this sequence is unmapped
SAM_MATEUNMAPPED  = 0x0008	# mate sequence is unmapped
SAM_QUERYSTRAND   = 0x0010	# query strand (0=>forward, 1=>reverse)
SAM_MATESTRAND    = 0x0020	# mate strand  (0=>forward, 1=>reverse)
SAM_PAIRFIRST     = 0x0040	# read is first in pair
SAM_PAIRSECOND    = 0x0080	# read is second in pair
SAM_NOTPRIMARY    = 0x0100	# alignment is not primary
SAM_CHECKFAIL     = 0x0200	# read failed platform/vendor checks
SAM_DUPLICATE     = 0x0400	# read is PCR or optical duplicate

class SAM_Reader(HTSeq.FileOrSequence):
    """
    Subclassed HTSeq.FileOrSequence which imitates what HTSeq.SAM_Reader
    does, except it yields our modified SAM_Alignment objects.
    """
    def __iter__(self):
        """
        Defined iterator method so that "object.next()" calls will work
        """
        for line in super(SAM_Reader,self).__iter__():
            if line.startswith('@'):
                continue
            algnt = SAM_Alignment(line)
            yield algnt

class SAM_Alignment(HTSeq.SAM_Alignment):
    """
    Subclass of the HTSeq.SAM_Alignment class.  This subclass stores
    additional information from the record such as the input record itself,
    and parses the tag elements which HTSeq.SAM_Alignment does not yet do.
    """

    def __init__(self, line):
        """
        Call super class __init__ and then supplement it by saving the
        input line and parsing the tags and saving them in a dict.
        """

        # Workarounds for HTSeq errors. 
        workaround_HTSeq_Exception_RNAME = False
        workaround_HTSeq_Exception_MRNM = False
        fields  = line.rstrip().split('\t')
        rname   = fields[2]
        mrnm    = fields[6]			
        flagint = int(fields[1])

        # Workaround for HTSeq "Malformed SAM line: RNAME != '*'..." ValueError 
        # exception in the super class __init__() method.  If we detect that
        # flag 0x0004 is set and the rname is not '*' we simply turn off that
        # flag before continuing, and then switch it back on once the super 
        # class __init__ has completed..   This may then result in another
        # ValueError exception because of a bad CIGAR string ('*') so we also 
        # change it and set it back after calling the super class __init__.

        if flagint & 0x0004 and rname != '*':
            workaround_HTSeq_Exception_RNAME = True
            flagint = flagint & ~0x0004
            fields[1] = '%d' % flagint
            cigar_string = fields[5]
            if fields[5] == '*':
                fields[5] = '%dM' % len(fields[9])
            line = '\t'.join(fields)

        # Workaround for HTSeq "Malformed SAM line: MRNM != '*'..." ValueError 
        # exception in the super class __init__() method.  In this case, SAM
        # files produced by novoalign v2.07.13 (previouly site was running 
        # v2.05.33) will not pass muster with HTSeq 0.4.4 if the MRNM is set
        # (i.e. is not '*') but 0x0008 flag is also set (MATE_UNMAPPED).
        # We account for that here by setting the MRNM to '*' manually
        if flagint & SAM_MATEUNMAPPED and mrnm != '*':
            workaround_HTSeq_Exception_MRNM = True
            fields[6] = '*'
            line = '\t'.join(fields)

        # Now call parent iitializer
        super(SAM_Alignment, self).__init__(line)

        # Undo changes made by RNAME workaround 
        if workaround_HTSeq_Exception_RNAME:
            flagint |= 0x0004
            # Save original flags integer
            self.flags = flagint
            fields[1] = '%d' % flagint
            fields[5] = cigar_string
            line = '\t'.join(fields)

        # Store parsed fields in an array we can avail of later
        self.fields = fields

        # Extract MRNM and tags fields directly from the input line
        lsplit = line.split("\t")
        self.tags = lsplit[6]
        self.__set_properties__()
        self.__set_additional_attributes__()


    def __getflagbit__(self, mask):
        """
        A convenient method of querying individual flag bits, this method 
        takes a single input argument: the bit (an index from 1 to 11).
          
        See also: __setflagbit__()
           
        Flag bits are defined as follows:
	SAM_NONE		0x0000	null
        SAM_MAPPED	 0	0x0001	read is paired in sequencing
        SAM_PROPERPAIR	 1	0x0002	read is mapped in proper pair
        SAM_UNMAPPED	 2	0x0004	this sequence is unmapped
        SAM_MATEUNMAPPED 3	0x0008	mate sequence is unmapped
        SAM_QUERYSTRAND	 4	0x0010	query strand (0=>forward, 1=>reverse)
        SAM_MATESTRAND	 5	0x0020	mate strand  (0=>forward, 1=>reverse)
        SAM_PAIRFIRST	 6	0x0040	read is first in pair
        SAM_PAIRSECOND	 7	0x0080	read is second in pair
        SAM_NOTPRIMARY		0x0100	alignment is not primary
        SAM_CHECKFAIL		0x0200	read failed platform/vendor checks
        SAM_DUPLICATE		0x0400	read is PCR or optical duplicate
        """
        return self.FLAG & mask is not 0
         

    def __setflagbit__(self, mask, state):
        """
        A convenient method of manipulating flag bits, this method takes
        two input arguments:  the bit (an index from 1 to 11) and the
        state 'on', 'off' or 'toggle'
           
        Flag bits are defined as follows:
	SAM_NONE		0x0000	null
        SAM_MAPPED	 0	0x0001	read is paired in sequencing
        SAM_PROPERPAIR	 1	0x0002	read is mapped in proper pair
        SAM_UNMAPPED	 2	0x0004	this sequence is unmapped
        SAM_MATEUNMAPPED 3	0x0008	mate sequence is unmapped
        SAM_QUERYSTRAND	 4	0x0010	query strand (0=>forward, 1=>reverse)
        SAM_MATESTRAND	 5	0x0020	mate strand  (0=>forward, 1=>reverse)
        SAM_PAIRFIRST	 6	0x0040	read is first in pair
        SAM_PAIRSECOND	 7	0x0080	read is second in pair
        SAM_NOTPRIMARY		0x0100	alignment is not primary
        SAM_CHECKFAIL		0x0200	read failed platform/vendor checks
        SAM_DUPLICATE		0x0400	read is PCR or optical duplicate
        """
        if state == 'on':
            self._flag |=  mask
        if state == 'off':
            self._flag &= ~mask
        if state == 'toggle':
            self._flag ^=  mask
         

    def __set_properties__(self):
        """
        This internal method sets up properties for each of the names fields
        in the SAM standard document.  Currently these are:

        QNAME
        FLAG
        RNAME
        POS
        MAPQ
        CIGAR
        MRNM
        MPOS
        ISIZE
        SEQ
        QUAL
        [TAG:VTYPE:VALUE,]

        We store the sequence and quality score string as read from the SAM 
        file, (HTSeq reverse-complements the sequence and phred scores for
        records which are deemed to have come from the reverse reference).
        """

        self._qname = self.read.name
        self._flag  = self.flags
        if hasattr(self.iv, 'chrom'):
            self._rname = self.iv.chrom
        else:
            self._rname = self.fields[2]
        self._pos = int(self.fields[3])
        self._mapq = self.aQual
        self._cigar = self.fields[5]
        self._mrnm = self.fields[6]
        self._mpos = int(self.fields[7])
        self._isize = self.inferred_insert_size
        # self.read is automaticaly reverse-complented if necessary. We use
        # self._read to store the read as found in the file (unchanged)
        self._read = self.read_as_aligned
        self._tags = []
        self._minPhred = None
        for tag in self.fields[11:]:
            self._tags.append(tag.split(':'))

    def __getMinPhred(self):
        """Compute minimum PHRED score and return it.  PHRED score is stored
        as a string of ASCII characters which encode the PHRED value plus 33
        (to make the characters fall within the normal printable ASCII range
        """
        val = inf
        for char in self.QUAL:
            val = min(ord(char)-33,val)
        return val

    def copy(self):
        return SAM_Alignment(self.STRING)

    def matchRNAME(self, rname):
        """
        Given an rname (string) returns True or False depending on whether 
        the rname (including wildcards) matches self.RNAME.  If rname is a tuple
        describing an rname string, minPos and maxPos, extract the first string
        and use that for matching
        """
        if type(rname) in (types.ListType, types.TupleType):
            if len(rname) > 0 and type(rname[0]) is types.StringType:
                rname = rname[0]
        if '*' in rname:
            rname = rname.replace('*', '.*')
            return re.search(rname, self.RNAME)
        else:
            return rname == self.RNAME

    def matchPOS(self, rangeTuple, flank, RHS=False):
        """
        Given a range tuple defined as (minPos,maxPos), returns True or False
        depending on whether self.POS falls inside the defined range, using
        "flank" as a margin.  If the optional RHS argument is True, match
        against the POSR2 position instead of the default POSL1 position.
        """
        minPos = rangeTuple[0] - flank
        maxPos = rangeTuple[1] + flank
        if minPos < 1: minPos = 1
        if maxPos == -1: maxPos = inf
        if RHS:
            return minPos <= self.POSR2 <= maxPos
        else:
            return minPos <= self.POSL1 <= maxPos

    def matchPHREDvalue(self, rangeTuple):
        """
        Given a range tuple defined as (minPhred,maxPhred), return True or 
        False depending on whether the phred score falls within the range
        """
        minPhred = rangeTuple[0]
        maxPhred = rangeTuple[1]
        if minPhred == -1: minPhred = inf
        if maxPhred == -1: maxPhred = inf
        return minPhred <= self.MINPHRED <= maxPhred

    def matchNMvalue(self, rangeTuple):
        """
        Given a range tuple defined as (minNM,maxNM), returns True or False
        depending on whether the NM:i tag integer value falls inside the 
        defined range.  If no NM:i tag exists, return False
        """
        minNM = rangeTuple[0]
        maxNM = rangeTuple[1]
        for t in self.TAGS: 
            if t[:2] == ['NM', 'i']:
                nmvalue = int(t[2])
                return minNM <= nmvalue <= maxNM
        return False

    def matchRPvalue(self, rangeTuple):
        """
        Given a range tuple defined as (minRP,maxRP), returns True or False
        depending on whether the ZN:i tag integer value falls inside the 
        defined range.  If no ZN:i tag or ZS:Z:R exists, assume a ZN:i value
        of zero and test for inclusion in minRP,maxRP range as usual (e.g.
        matchRPvalue will pass (return True) if no ZS:Z:R/ZN:i tags are present
        but minRP is 0)
        """
        RPvalue = 0
        minRP = rangeTuple[0]
        maxRP = rangeTuple[1]
        for t in self.TAGS: 
            if t[:2] == ['ZN', 'i']:
                RPvalue = int(t[2])
                break
        return minRP <= RPvalue <= maxRP

    def matchRNAMEandPOS(self, rnameTuple, flank, RHS=False):
        """
        Given an RNAME tuple defined as (RNAME, minPos, maxPos), returns True 
        or False depending on whether self.RNAME matches and self.POS falls 
        inside the defined range.  The optionsl RHS will be passed to matchNAME,
        where, if True, the POSR2 value will be used for matching instead of the
        default posL1
        """
        (rname, minPos, maxPos) = rnameTuple
        return self.matchRNAME(rname) and self.matchPOS((minPos,maxPos), 
                                                        flank, RHS)

    @property
    def QNAME(self):
        """
        SAM record QNAME value
        """
        return self._qname
        
    @property
    def FLAG(self):
        """
        SAM record FLAG value
        """
        return self._flag

    @property
    def RNAME(self):
        """
        SAM record RNAME value
        """
        return self._rname

    @property
    def POS(self):
        """
        SAM record POS value
        """
        return self._pos

    @property
    def POSL1(self):
        """
        L1 (left-hand) position POS value (always same as POS)
        """
        return self.POS

    @property
    def POSR2(self):
        """
        R1 (left-hand) position POS value (always POS + LENGTH)
        (See samkit 4 spec 4)
        """
        return self.POS + self.LENGTH

    @property
    def READ_BEGIN(self):
        """
        Overall position of first base in the sequence. Same as POSL1
        """
        ### NOTE: This was implemented in SBLsamkitEFE 4.15.  Not required since
        ### the definition of POSR2 has been corrected??
        return self.POS

    @property
    def READ_END(self):
        """
        Overall position of last base in the sequence. BEGIN + readLength
        """
        ### NOTE: This was implemented in SBLsamkitEFE 4.15.  Not required since
        ### the definition of POSR2 has been corrected??
        return self.READ_BEGIN + self.LENGTH

    @property
    def LENGTH(self):
        """
        Sequence length
        """
        return len(self.SEQ)

    @property
    def MAPQ(self):
        """
        SAM record MAPQ value
        """
        return self._mapq

    @property
    def CIGAR(self):
        """
        SAM record CIGAR value
        """
        return self._cigar

    @property
    def MRNM(self):
        """
        SAM record MRNM value
        """
        return self._mrnm

    @property
    def MPOS(self):
        """
        SAM record MPOS value
        """
        return self._mpos

    @property
    def ISIZE(self):
        """
        SAM record ISIZE value
        """
        return self._isize

    @property
    def SEQ(self):
        """
        SAM record unmoodified SEQ string (*NEVER* reverse complemented)
        """
        return self._read.seq

    @property
    def QUAL(self):
        """
        SAM record unmodified QUAL string
        """
        return self._read.qualstr

    @property
    def MINPHRED(self):
        """
        Minimim Phred score
        """
        if self._minPhred is None:
            self._minPhred = self.__getMinPhred()
        return self._minPhred

    @property
    def TAGS(self):
        """
        SAM record tags array:
        """
        return self._tags

    @property
    def STRAND(self):
        """
        SAM record strand
        """
        return self._strand

    @property
    def mapped(self):
        """
        Sam record is marked as mapped.  If the CIGAR string is '*' we can 
        assume the record is not uniquely aligned.  Don't confuse this with
        the "mapped" flag bits.
        """
        return self.CIGAR != '*'

    @property
    def STRING(self):
        """
        This property returns a string format of the current fields in the
        SAM_Alignment object.  So at any time the state of the SAM record
        (e.g. after altering attributes e.g. manual reverse-complement) can be
        converted back into a string.  By default self.SEQ will return the
        sequence as read in.
        """
        samrecord = ['%s' % self.QNAME,
                     '%d' % self.FLAG,
                     '%s' % self.RNAME,
                     '%d' % self.POS,
                     '%d' % self.MAPQ,
                     '%s' % self.CIGAR,
                     '%s' % self.MRNM,
                     '%d' % self.MPOS,
                     '%d' % self.ISIZE,
                     '%s' % self.SEQ,
                     '%s' % self.QUAL,
                        ]
        for tag in self.TAGS:
            samrecord.append(':'.join(tag))
        return '\t'.join(samrecord)

    #---------------- L: Record location mapping properties
    @property
    def M(self):
        """
        Record is mapped - either uniquely mapped or repeat mapped
        """
        return self.U or self.rM

    @property
    def U(self):
        """
        A convenience property which returns True if this record is uniquely
        mapped
        """
        return self.aligned and self.mapped and not self.rM

    @property
    def nM(self):
        """
        A convenience property which returns True if this record is not mapped
        """
        if self.U or self.M:
            return False

        # Novoalign flags indicate record which isn't mapped
        if ['ZS', 'Z', 'NM'] in self.TAGS:
            return True

        # Per Nick's description:
        #
        # stampy output flags won't have ['ZS', 'Z', 'NM'] in the tags, but
        # will have one of the following flag values:
        # 
        #    77/141 flags to detect pairs where both reads are unmatched
        #    73/133 (UnM PE pair where second read in pair is unmapped and 
        #            mapped read is in forward orientation) 
        #    89/165 (UnM PE pair where second read in pair is unmapped and 
        #            mapped read is in reverse orientation) 
        #    69/137 (UnM PE pair where first read in pair is unmapped and 
        #            mapped read is in forward orientation) 
        #   101/153 (UnM PE pair where first read in pair is unmapped and 
        #            mapped read is in reverse orientation) 
        #
        # Instead of just comparing numerical values we'll test the mask bit,
        # which shows as expected the following flag values indicate an unmapped
        # sequence:
        #    77, 141, 133, 165, 69, 101
        return self.__getflagbit__(SAM_UNMAPPED)

    @property
    def rN(self):
        """
        A convenience property which returns True if this record is 
        repeat mapped (Novoalign)
        """
        return ['ZS', 'Z', 'R'] in self.TAGS
        
    @property
    def rS(self):
        """
        A convenience property which returns True if this record is 
        repeat mapped (stampy)
        """
        result = False
        for tag in self.TAGS:
            if len(tag) > 2 and tag[:2] == ['XA', 'Z']:
                result = True
                break
        return result

    @property
    def rM(self):
        """
        Top-level repeat-mapped filter (matches rN and rS)
        """
        return self.rN or self.rS

    @property
    def rMM(self):
        """
        Convenience property - used by matchWithRangeList in Stream.py. Returns
        true if the record is marked repeat-mapped, but *is* mapped (has an
        RNAME and POS
        """
        return self.aligned and self.mapped and self.rM

    @property
    def X(self):
        """
        A convenience property - returns True if *either* not mapped or repeat
        mapped
        """
        return not self.U

    @property
    def A(self):
        """
        A convenience property which always returns True (match Any)
        mapped
        """
        return True
    
    @property
    def QC(self):
        """
        A convenience property which returns True if this record has been 
        marked quality failure by the aligner ( 'ZS:Z:QC' )
        """
        return ['ZS', 'Z', 'QC'] in self.TAGS

    @property
    def CLC1(self):
        """
        This property returns True if self.CIGAR contains the letter 'I' at the
        beginning, after up to three numerical digits, or at the end.  See the
        --ssN,CLC1 filter for more information
        """
        if re.match('^[0-9]{0,3}I', self.CIGAR) is not None or \
           re.match('.*I$', self.CIGAR) is not None:
            return True
        return False
 
    #---------------- s: Record strand mapping properties
    @property
    def f(self):
        """
        Convenience property which returns True if the record is forward mapped
        """
        return self.M and self.STRAND == '+'

    @property
    def r(self):
        """
        Convenience property which returns True if the record is reverse mapped
        """
        return self.M and self.STRAND == '-'

    @property
    def x(self):
        """
        Returns True if record is not uniquely mapped
        """
        return not self.U

    @property
    def s(self):
        """
        Returns True is either forward or reverse mapped
        """
        return self.f or self.r


    def __set_additional_attributes__(self):
        """
        Sets up additional useful attributes, outside the standard SAM spec.
        These are used as convenience by samkit3.
        """
        self.__set_strand__()

    def __set_strand__(self):
        """
        Sets the strand direction where applicable.  If the strand is mapped
        (0x0004 is not set, 0x0008 is), use flag 0x0010.  If the strand is 
        unmapped but the mate is mapped, (0x0004 is set, 0x0008 is not), 
        instead use flag 0x00020
	SAM_NONE		0x0000	null
        SAM_MAPPED		0x0001	read is paired in sequencing
        SAM_PROPERPAIR		0x0002	read is mapped in proper pair
        SAM_UNMAPPED		0x0004	this sequence is unmapped
        SAM_MATEUNMAPPED	0x0008	mate sequence is unmapped
        SAM_QUERYSTRAND		0x0010	query strand (0=>forward, 1=>reverse)
        SAM_MATESTRAND		0x0020	mate strand  (0=>forward, 1=>reverse)
        SAM_PAIRFIRST		0x0040	read is first in pair
        SAM_PAIRSECOND		0x0080	read is second in pair
        SAM_NOTPRIMARY		0x0100	alignment is not primary
        SAM_CHECKFAIL		0x0200	read failed platform/vendor checks
        SAM_DUPLICATE		0x0400	read is PCR or optical duplicate
        """
        if not self.__getflagbit__(SAM_UNMAPPED):
            if self.__getflagbit__(SAM_QUERYSTRAND):
                self._strand = '-'
            else:
                self._strand = '+'
        elif self.__getflagbit__(SAM_UNMAPPED) and \
             not self.__getflagbit__(SAM_MATEUNMAPPED):
            if self.__getflagbit__(SAM_MATESTRAND):
                self._strand = '-'
            else:
                self._strand = '+'
        else:
            self._strand = '?'


def pair_SAM_alignments(alignments, verbose=False ):
    """
    Basically a duplication of HTSeq.pair_SAM_alignments() function which
    suppresses some warnings, and utilizes our SAM_Reader and SAM_Alignment
    classes.  Adds a second optional argument 'verbose' which by default is
    False, to quence verbose warning output
    """

    def check_is_pe( read ):
        """
        Identical to the check_is_pe function embdedded in the original 
        HTSeq.pair_SAM_alignments() function.
        """
        if not read.paired_end:
            raise ValueError("%s is not a paired_end read" % repr(read))

    def process_paired_reads(read1,read2):   
        """
        An enhanced implementation of the process_paired_reads() function
        embedded in the original HTSeq.pair_SAM_alignments() function.  This
        implementation removes some spurious warning messages, and also
        checks for the existance of certain object attributes such as 
        iv.start_as_pos, before comparing them.
        """
        if read1.pe_which == "second":
            aux = read1
            read1 = read2
            read2 = aux
        if not ((read1.pe_which == "first" and read2.pe_which == "second" ) or 
                (read1.pe_which == "unknown" and read2.pe_which == "unknown")):
            if verbose:
                warnings.warn("Incorrect first/second assignments in mate" +
                              " pairs " + read1.read.name)
        if not ( read1.proper_pair and read2.proper_pair ):
            if verbose:
                warnings.warn("Incorrect 'proper_pair' flag value for read" +
                              " pair " + read1.read.name)
        if hasattr(read2.iv, "start_as_pos"):
            if not (read1.mate_start == read2.iv.start_as_pos and 
                    read2.mate_start == read1.iv.start_as_pos ):
                if verbose:
                    warnings.warn("Read pair " + read1.read.name +
                                  " show inconsistency between 'iv'" +
                                  " and 'mate_start' values" )
        return (read1, read2)

    def process_single_read( read ):
        if read.mate_aligned:         
            warnings.warn("Read " + read.read.name + 
                          " claims to have an aligned mate which could not" +
                          " be found. (Is the SAM file properly sorted?)")
        if read.pe_which == "second":
            return (None, read)
        else:
            return (read, None)

    def names_match(read1, read2):
        """Compare name of both records.  If it's a match, return True, all is good.
        If not however, compare the two strings again after stripping off any trailing
        _1 or _2 or /1 or /2 strings (from fastqkit)."""

        if read1.read.name == read2.read.name:
            return True
        name1 = re.sub('[/_][12]$', '', read1.read.name)
        name2 = re.sub('[/_][12]$', '', read2.read.name)
        if name2 == name2:
            return True
        return False

    # Begin operation of main method #
    alignments = iter(alignments)
    read1 = None
    read2 = None
    while True:
        if read1 is None:
            read1 = alignments.next()
            check_is_pe(read1)
        read2 = alignments.next()
        check_is_pe(read2)      
        if names_match(read1, read2):
            yield process_paired_reads(read1, read2)
            read1 = None
            read2 = None
        else:
            yield process_single_read(read1)
            read1 = read2
            read2 = None

