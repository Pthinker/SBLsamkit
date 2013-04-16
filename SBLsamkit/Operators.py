####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""Constants defining operators and location and strand selector strings"""

#============================================================================#
# IMPORTANT.  Please refer to the SAMStream.evaluatePairMapping() object method
# defined in Stream.py.  If you add or remove and Location or strand mapping
# names in this file (LL_ or SS_), YOU MUST UPDATE the aforementioned method
# with the same changed names.
#============================================================================#

# Location mapping names. See LL properties section of SAM.SAMPair class
LL_UX      = 'UX'
LL_UNM     = 'UnM'
LL_URM     = 'UrM'
LL_URN     = 'UrN'
LL_URS     = 'UrS'
LL_UU      = 'UU'
LL_UUE     = 'UUE'
LL_UUD     = 'UUD'
LL_UA      = 'UA'
LL_XX      = 'XX'
LL_XA      = 'XA'
LL_MX      = 'MX'
LL_MM      = 'MM'
LL_MNM     = 'MnM'
LL_MRM     = 'MrM'
LL_MRN     = 'MrN'
LL_MRS     = 'MrS'
LL_MME     = 'MME'
LL_MMD     = 'MMD'
LL_MA      = 'MA'
LL_NMNM    = 'nMnM'
LL_RMRM    = 'rMrM'
LL_RNRN    = 'rNrN'
LL_RSRS    = 'rSrS'
LL_NMRM    = 'nMrM'
LL_NMRN    = 'nMrN'
LL_NMRS    = 'nMrS'
LL_NMA     = 'nMA'
LL_RMA     = 'rMA'
LL_RNA     = 'rNA'
LL_RSA     = 'rSA'
LL_QC      = 'QC'
LL_CLC1    = 'CLC1'
LL_AA      = 'AA'

# Strand mapping names. See SS properties section of SAM.SAMPair class
SS_FR      = 'fr'
SS_RF      = 'rf'
SS_SSD     = 'ssd'
SS_FX      = 'fx'
SS_RX      = 'rx'
SS_FF      = 'ff'
SS_RR      = 'rr'
SS_SSE     = 'sse'
SS_AA      = 'aa'

# Filter names
FP_NMF     = 'NMf'
FP_NMFA    = 'NMfa'
FP_RPF     = 'RPf'
FP_RPFA    = 'RPfa'
FP_PHREDF  = 'phredf'
FP_PHREDFA = 'phredfa'

# Operator names
OP_FASTQ   = 'fastq'
OP_SAM     = 'sam'
OP_SAMPP   = 'sam,pp'
OP_FASTQPP = 'fastq,pp'
OP_PP      = 'pp'
OP_SH      = 'sh'
OP_ORDER   = 'order'
OP_SE      = 'se'
OP_NOUT    = 'nout'
OP_INFO    = 'info'
OP_NOINFO  = 'noinfo'
OP_GZ      = 'gz'
OP_CSV     = 'csv'
OP_PASS    = 'pass'
OP_EX      = 'ex'
OP_SAMPLE  = 'sample'

# Order secondary operators (order mode)
ORD_LL     = 'LL'
ORD_SS     = 'ss'
ORD_POS    = 'pos'
ORD_RNG    = 'rng'
ORD_RNGP   = 'rngp'

# Range operators and subcommands
RNG_RNG        = 'rng'
RNG_RNGF       = 'rngf'
RNG_RNGPAIR    = 'rngp'
RNG_RNGPAIRF   = 'rngpf'
RNG_RANGES     = (RNG_RNG, RNG_RNGPAIR)
RNG_RANGEFILES = (RNG_RNGF, RNG_RNGPAIRF)
RNG_FLANK      = 'flank'
RNG_FOOTPRINT  = 'footprint'
RNG_REPEAT     = 'repeats'
RNG_COMMANDS   = (RNG_FLANK, RNG_FOOTPRINT, RNG_REPEAT)
RNG_FP_L1      = 'L1'
RNG_FP_L1ANDR2 = 'L1ANDR2'
RNG_FP_L1ORR2  = 'L1ORR2'
RNG_FP_R2      = 'R2'
RNG_FOOTPRINTS = (RNG_FP_L1, RNG_FP_L1ANDR2, RNG_FP_L1ORR2, RNG_FP_R2)

# Recognized streams
STREAMS    = (('st', 
             'suffix,[LL,][ss,][NMf,min,max,][RPf,min,max,][phredf,min,max,][iL,min,max,][[rng (rangelist),|rngp (range pair list),|rngf,rangefilename,|rngpf,rangepairfilename,[flank,value,][footprint,value,][repeats]][,operator1[,operator2...]]'),
             )

# Generates a list of Location mapping names from LL_*
LLnames    = [ LL_UX,
               LL_UNM,
               LL_URM,
               LL_URN,
               LL_URS,
               LL_UU,
               LL_UUE,
               LL_UUD,
               LL_UA,
               LL_XX,
               LL_XA,
               LL_MX,
               LL_MM,
               LL_MNM,
               LL_MRM,
               LL_MRN,
               LL_MRS,
               LL_MME,
               LL_MMD,
               LL_MA,
               LL_NMNM,
               LL_RMRM,
               LL_RNRN,
               LL_RSRS,
               LL_NMRM,
               LL_NMRN,
               LL_NMRS,
               LL_NMA,
               LL_RMA,
               LL_RNA,
               LL_RSA,
               LL_QC,
               LL_CLC1,
               LL_AA,
             ]

# Generates a list of strand mapping names from SS_*
SSnames    = [ SS_FR,
               SS_RF,
               SS_SSD,
               SS_FX,
               SS_RX,
               SS_FF,
               SS_RR,
               SS_SSE,
               SS_AA,
             ]

# Generates a list of filter names from FP_*
FPnames    = [ FP_NMF,
               FP_NMFA,
               FP_RPF,
               FP_RPFA,
               FP_PHREDF,
               FP_PHREDFA,
             ]

# Generates a list of operator names from OP_*
OPnames    = [ OP_FASTQ,
               OP_SAM,
               OP_FASTQPP,
               OP_SAMPP,
               OP_PP,
               OP_SH,
               OP_ORDER,
               OP_SE,
               OP_NOUT,
               OP_INFO,
               OP_NOINFO,
               OP_GZ,
               OP_CSV,
               OP_PASS,
               OP_EX,
               OP_SAMPLE,
             ]

# Generates a list of order operator modes from ORD_*
ORDnames   = [ ORD_LL,
               ORD_SS,
               ORD_POS,
               ORD_RNG,
               ORD_RNGP,
             ]

# Operator help text
OPERATORHELP = ( 
  'Sub-stream specifiers:',
  '',
  '  %s' % OP_FASTQ,
  '  %s' % OP_SAM,
  '',
  'Operators - these act on the stream to their left:',
  '',
  '  %s' % OP_SH,
  '  %s' % OP_ORDER,
  '  %s' % OP_SE,
  '',
  '  %-12s %s' % (OP_NOUT,  'no file(s) output'),
  '  %-12s %s' % (OP_INFO,  'like nout, but info file(s) are still written'),
  '  %-12s %s' % (OP_NOINFO,'no info file(s)'),
  '  %-12s %s' % (OP_GZ ,   'file(s) output as compressed file(s)'),
  '  %-12s %s' % ('%s,N,M' % OP_SAMPLE,
                            'only accept first N of each M matches'),
  '  %-12s %s' % ('%s,N' % OP_CSV, 
                            'create or add range qualifying records to .N.csv'),
  '',
  'Modifiers:',
  '',
  '  %s\t%s' % (OP_PASS,'causes all records to be passed onto the next stage'),
  '  %s\t%s' % (OP_EX, 
            'parsed records are passed to the next stage, others are streamed'),
  '',
  'L = record location type:',
  '',
  '  M     :mapped - has an RNAME and a pos',
  '  U     :uniquely mapped - RNAME and pos present, no ZS:Z:R field',
  '  nM    :not mapped - no RNAME and pos, no  ZS:Z:R field',
  '  rM    :repeat mapped. All records with ZS:Z:R, even if no RNAME & pos',
  '  X     :not U',
  '  A     :any',
  '',
  'LL = record location type of a PE read. note that some values are supersets',
  '     of the others.',
  '',
  '  %s' % LL_UX,
  '  %s' % LL_UNM,
  '  %s' % LL_URM,
  '  %s' % LL_UU,
  '  %s\t%s' % (LL_UUE, ':both U, mapping to the same RNAME'),
  ' ^%s\t%s' % (LL_UUE, ':not UUE'),
  '  %s\t%s' % (LL_UUD, ':both U, mapping to different RNAMEs'),
  '  %s\t%s' % (LL_UA,  ':at least one U'),
  '  %s' % LL_MX,
  '  %s' % LL_MNM,
  '  %s' % LL_MRM,
  '  %s' % LL_MM,
  '  %s' % LL_MME,
  ' ^%s\t%s' % (LL_MME, ':not MME'),
  '  %s\t%s' % (LL_MMD, ':both M, mapping to different RNAMEs'),
  '  %s\t%s' % (LL_MA,  ':at least one record M'),
  '  %s' % LL_NMNM,
  '  %s' % LL_RMRM,
  '  %s' % LL_NMRM,
  '  %s' % LL_NMA,
  '  %s' % LL_RMA,
  '',
  '  %s\t%s' % (LL_XX, ':at least on record is X'),
  '  %s\t%s' % (LL_AA, ':all records'),
  '  %s\t%s' % (LL_QC, 
          ':a record pair which has failed any QC test (not yet implemented)'),
  '  %s\t%s' % (LL_CLC1,
          ': a record pair which is incompatible with SAM input into CLC Bio.'),
  '',
  'ss = strand mapping type of a PE read:  second column below is default LL',
  '     field where applicable:',
  '',
  '  s types:',
  '  f	:record maps to forward strand',
  '  r	:record maps to reverse strand',
  '  x	:record not uniquely mapped',
  '  s	:record is either f or r',
  '',
  '  ss:\tdefault LL',
  '  %s\t%s' % (SS_FR,  ':MME'),
  ' ^%s\t%s' % (SS_FR,  ':AA      :not fr'),
  '  %s\t%s' % (SS_RF,  ':MME'),
  '  %s\t%s' % (SS_SSD, ':MME     :read is either fr or rf'),
  '  %s\t%s' % (SS_FX,  ':MX'),
  '  %s\t%s' % (SS_RX,  ':MX'),
  '  %s\t%s' % (SS_FF,  ':MM'),
  '  %s\t%s' % (SS_RR,  ':MM'),
  '  %s\t%s' % (SS_SSE, ':MM      :record is either ff or rr'),
  '  %s\t%s' % (SS_AA,  ':AA      :all records'),
  '')
