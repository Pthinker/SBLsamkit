####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####



1) INSTALLATION:

To install the SBLsamkit module and commands just issue the following command.  
You will likely need to execute this command as root, or using sudo:

    python setup.py install

The SBLsamkit module files will be installed in the python "site-packages" 
directory while the executable commands (see USAGE, below) are copied into the 
same directory as the python executable (typically /usr/bin or /usr/local/bin).

If installing on Windows however, the executable commands are installed in the 
Python "Scripts" subdirectory.  You may need to add the "Scripts" subdirectory 
to your Path environment -- see "Using Python on Windows" at 
http://docs.python.org/using/windows.html for instructions on how to do that.



2) USAGE:

The SBLsamkit package currently ships with 1 command:

Usage: samkit.py [options]

Options:
  --version          show program's version number and exit
  -h, --help         show this help message and exit
  --if=INPUTFILE     input filename                  [REQUIRED]
  --limit=LIMIT      limit to N record pairs         (default: no limit)
  --dr=OUTPUT_DIR    Specify output directory        (default: automatic)
  --np=NPROC         number of parallel processes    (default: 1 (single))
  --chunk=CHUNKSIZE  number of record pairs per chunk (default: 50)
  --unsorted         Disable output record sorting when using multiple
                     processors (default: don't)
  --hwm=QHIGH        High watermark num objects in queue (default:500000)
  --lwm=QLOW         Low watermark num objects in queue (default:70% QHIGH)
  --ini=INI          Read predefined set of filters from INI file(s)
  --stats=STATS      Select [f]ull or [c]ompact .info output
  --percentile=PCT   Set ranges for percentile splits

Sub-stream specifiers:

  fastq
  sam

Operators - these act on the stream to their left:

  sh
  order
  se

  nout         no file(s) output
  info         like nout, but info file(s) are still written
  noinfo       no info file(s)
  gz           file(s) output as compressed file(s)
  sample,N,M   only accept first N of each M matches
  csv,N        create or add range qualifying records to .N.csv

Modifiers:

  pass	causes all records to be passed onto the next stage
  ex	parsed records are passed to the next stage, others are streamed

L = record location type:

  M     :mapped - has an RNAME and a pos
  U     :uniquely mapped - RNAME and pos present, no ZS:Z:R field
  nM    :not mapped - no RNAME and pos, no  ZS:Z:R field
  rM    :repeat mapped. All records with ZS:Z:R, even if no RNAME & pos
  X     :not U
  A     :any

LL = record location type of a PE read. note that some values are supersets
     of the others.

  UX
  UnM
  UrM
  UU
  UUE	:both U, mapping to the same RNAME
 ^UUE	:not UUE
  UUD	:both U, mapping to different RNAMEs
  UA	:at least one U
  MX
  MnM
  MrM
  MM
  MME
 ^MME	:not MME
  MMD	:both M, mapping to different RNAMEs
  MA	:at least one record M
  nMnM
  rMrM
  nMrM
  nMA
  rMA

  XX	:at least on record is X
  AA	:all records
  QC	:a record pair which has failed any QC test (not yet implemented)
  CLC1	: a record pair which is incompatible with SAM input into CLC Bio.

ss = strand mapping type of a PE read:  second column below is default LL
     field where applicable:

  s types:
  f	:record maps to forward strand
  r	:record maps to reverse strand
  x	:record not uniquely mapped
  s	:record is either f or r

  ss:	default LL
  fr	:MME
 ^fr	:AA      :not fr
  rf	:MME
  ssd	:MME     :read is either fr or rf
  fx	:MX
  rx	:MX
  ff	:MM
  rr	:MM
  sse	:MM      :record is either ff or rr
  aa	:AA      :all records


Streams, N in range 1-99: 
  --stN=              suffix,[LL,][ss,][NMf,min,max,][RPf,min,max,][phredf,min,m
                    ax,][iL,min,max,][[rng (rangelist),|rngp (range pair list)
                    ,|rngf,rangefilename,|rngpf,rangepairfilename,[flank,value
                    ,][footprint,value,][repeats]][,operator1[,operator2...]]




3) RUNTIME CONFIGURATION

# <------- SAMPLE 'samkit.ini' configuration file ----------- [cut here] ------>
#
# Beginning on this line is an example configuration file for samkit.  It can 
# be used to describe an ordered sequence of stream filters, as shown below.
# Cut-and-paste the sections below into a new file which must reside in the
# current working directory when running samkit.py.  Then you can call 
# samkit.py without specifying streams in the normal way, but instead by 
# calling one or more sections of this file by name with the --ini argument. 
#
# Example:
#
#   samkit.py --if=xmr002.sam --ini=samkit.ini,test1
#
# The implied command line arguments (i.e. as determined by the options in a 
# given section) will be recorded in any .info file which is output.  You may
# pass the --ini argument multiple times, and each time request a different
# section of the ini file, and note that you can also use multiple ini files.
# Be careful however that no two sections contain a stream of the same name,
# otherwise an error will be reported.
#
# The section name may be omitted, in which case the standard section called
# 'section01' will be used by default.
#
# Example:
#
#   samkit.py --if=xmr002.sam --ini=samkit1.ini,test1 --ini=samkit2.ini
#
# This will configure all the streams in section 'test1' of the samkit1.ini file
# and section 'section01' of the samkit2.ini file.
#
# More detailed examples of samkit.ini files are provided in the ./examples
# directory.

[test1]
st1=UUE,fr,stp

[test2]
st2=UU,ff,0,-1,rjp.ff,sam
st3=UU,rr,0,-1,rjp.rr,sam

#
# <------- SAMPLE 'samkit.ini' configuration file ----------- [cut here] ------>


4) API DOCUMENTATION

The SBLsamkit module is documented throughout with in-code comments.  These 
have been parsed by the "epydoc" utility (epydoc.sourceforge.net) to generate a 
comprehensive (hopefully!) set of API dpcumentation for developers.  You'll find
this in the docs/developer directory. Start with the "docs/developer/index.html"
file.

