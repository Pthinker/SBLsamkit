#!/usr/bin/env python
####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

import sys
from SBLsamkit import __version__
from SBLsamkit.Operators import *
from SBLsamkit.Parser import OptionParser, parse_stream, parse_int, parse_nproc
from SBLsamkit.Processor import ProcessStreams
from SBLsamkit.Stream import ConfigureStreams
from SBLsamkit.Utils import Error, ErrorHandler

version = "%prog " + "%s" % __version__

# Configure for unbuffered stdin,stdout,stderr (i.e. "python -u")
PYTHONUNBUFFERED = 'yes'

# Set up option parser and add each option defined in the options dictionary,
# in the order they are stored in the optionkeys list
parser = OptionParser(version=version)

# Add global command line options.  Note that for some integer types we 
# set the type as 'string' and use the parse_int callback.  This callback
# will accept ints with multipliers such as 100K or 25M (case-insensitive) and
# will store the resulting value as an int
parser.add_option('--if',
                  action='store',
                  dest='inputfile',
                  type='string',
                  default=None,
                  help="input filename                  [REQUIRED]")
parser.add_option('--limit', 
                  action='callback',
                  dest='limit',
                  callback=parse_int,
                  type='str',
                  default=None,
                  help="limit to N record pairs         (default: no limit)")
parser.add_option('--dr',
                  action='store',
                  dest='output_dir',
                  type='str',
                  default=None,
                  help='Specify output directory        (default: automatic)')
parser.add_option('--np',
                  action='callback',
                  dest='nproc',
                  callback=parse_nproc,
                  type='str',
                  default=1,
                  help='number of parallel processes    (default: 1 (single))')
parser.add_option('--chunk',
                  action='store',
                  dest='chunksize',
                  type='int',
                  default=50,
                  help='number of record pairs per chunk (default: 50)')
parser.add_option("--unsorted",
                  action="store_true",
                  default=False,
                  dest="mp_unsorted",
                  help="Disable output record sorting when using multiple processors (default: don't)")
parser.add_option('--hwm',
                  action='store',
                  dest='QHIGH',
                  type='int',
                  default=500000,
                  help='High watermark num objects in queue (default:500000)')
parser.add_option('--lwm',
                  action='store',
                  dest='QLOW',
                  type='int',
                  default=350000,
                  help='Low watermark num objects in queue (default:70% QHIGH)')
parser.add_option('--ini',
                  action='append',
                  dest='ini',
                  type='str',
                  default=None,
                  help='Read predefined set of filters from INI file(s)')
parser.add_option('--stats', 
                  action='store',
                  dest='stats',
                  default='full',
                  help='Select [f]ull or [c]ompact .info output')
parser.add_option('--percentile',
                  action='store',
                  dest='pct',
                  type='str',
                  default=None,
                  help='Set ranges for percentile splits')

# Add supported streams
for sname, sdoc in STREAMS:
    parser.add_multi_option(sname,
                            action='callback',
                            callback=parse_stream,
                            dest=sname,
                            type='string',
                            help=sdoc)

# If we have NO arguments or options, default to displaying help output:
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(-1)

# Parse input arguments use the older exception syntax to support Python 2.5.
try:
    opts,args,order,extra,inifiles = parser.parse_args()
    #if not hasattr(opts, 'nproc'):
    #    opts.nproc = 1
except:
    ErrorHandler()

# Configure streams we've just set up based on certain input options. This
# will also raise an exception is a stream has been declared more than once
ConfigureStreams(opts)

# Make sure no unrecognized arguments remain
if len(args):
    Error("Unrecognized arguments: %s" % repr(args))

# Process streams and catch Ctrl-C if pressed
try:
    ProcessStreams(opts, order, __file__, extra, inifiles)
except KeyboardInterrupt:
    print "Process interrupted by user"
except:
    raise

