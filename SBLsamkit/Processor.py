####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""Main driver functions for samkit4"""

import os
import re
import csv
import sys
import time
import gzip
import datetime
from Info import InfoFile
from SAM import ReadSAMFile
from Operators import *
from MultiProcess import *
from Stream import SAMStream
from Stats import *
from Utils import *


def validateOptions(opts):
    """
    Validates input arguments are correct, for example making sure the
    requested input file exists, limit arguments are correct, and so on
    """

    # Translate tilde at head of input filename and output dir if applicable
    opts.inputfile  = ExpandPath(opts.inputfile)
    opts.output_dir = ExpandPath(opts.output_dir)

    # Set up InfoFile object
    opts.infofile   = InfoFile(filename=opts.inputfile)

    # Verify input file has been provided and is available
    if opts.inputfile is None:
        Error("You must provide an input file (--if)")
    if not FileExists(opts.inputfile):
        Error("Input file %s not found" % opts.inputfile)
    if not FileReadable(opts.inputfile):
        Error("Cannot read input file %s" % opts.inputfile)

    # Verify stats argument is valid
    if not opts.stats.startswith('f') and not opts.stats.startswith('c'):
        Error("Invalid stats mode: %s" % opts.stats)
    else:
        opts.stats = opts.stats[0].lower()

    # Process csv,N operators
    opts.csv = []
    for stream in opts.orderedStreams:
        if stream.op(OP_CSV): 
            opts.csv.append({'stream' : stream.name,
                             'suffix' : stream.ops['csvExt'],
                             'op'     : OP_CSV,
                             'rnames' : stream.ops['rnames']})

    # Process sample,N,M operators
    for stream in opts.orderedStreams:
        if stream.op(OP_SAMPLE): 
            opts.csv.append({'stream' : stream.name,
                             'suffix' : stream.ops['csvExt'],
                             'op'     : OP_CSV,
                             'rnames' : stream.ops['rnames']})


def makeDirectories(opts, callername):
    """
    Create output directories
    """

    now = datetime.datetime.now().strftime('%Y.%m.%d_%H.%M.%S')
    opts.outputDirectories = []

    outputdir = opts.output_dir
    if outputdir is None:
        caller = os.path.splitext(os.path.split(callername)[1])[0]
        outputdir = now + '_' + caller

    opts.outputDirectories = [outputdir,]

    for directory in opts.outputDirectories:
        os.mkdir(directory)


def cleanDirectories(opts):
    """
    Remove empty output directories
    """
    for directory in opts.outputDirectories:
        pass
        try:
            os.rmdir(directory)
        except OSError:
            pass


def reportStats(opts, total, heartbeat):
    """
    Compute all non-stream related statistics for a completed run
    and store the output in the InfoFile object (referenced in 
    opts.infofile).  This information is shared among all .info files. Stats
    specific to a given stream or run can be set in the writeInfoFiles
    function as each stream is iterated over.
    """
    opts.infofile.stats['runtime']  = heartbeat.TotalRunTimeInSeconds
    opts.infofile.stats['count']    = heartbeat.count

    errors   = []
    warnings = []
    for stream in opts.orderedStreams:
        errors.extend(stream.errors)
        warnings.extend(stream.warnings)
    opts.infofile.stats['errors']   = errors
    opts.infofile.stats['warnings'] = warnings

    # The "generalStats" gets passed directly to the InfoFile append2ELN()
    # procedure, and the streamSTATS gets returned.
    # where each element is a 3-item tuple of values,  e.g:
    #   (".uue", "st1", 327)
    generalSTATS, streamSTATS = report_general_stats(opts, total)
    opts.infofile.append2ELN(generalSTATS)
    opts.infofile.append2ELN(report_stream_summary(opts))
    opts.infofile.append2ELN(report_rname_summary(opts))
    opts.infofile.append2ELN(report_rname_stats(opts))
    if opts.stats == 'f':
        opts.infofile.append2ELN(report_size_distributions(opts))
    return streamSTATS


def writeInfoFiles(opts, outputdir, finalOutput, streamSTATS, extraFiles=[]):
    """
    Generate .info file(s).   
    """

    if not outputdir.endswith('/'):
        outputdir += '/'

    AllFilesWritten = []
    AllRangeFiles   = []
    for stream in opts.orderedStreams:
        for filename in stream.OutputFiles:
            AllFilesWritten.append(filename)
        for rngfile in stream.ops['rngfile']:
            if rngfile not in AllRangeFiles:
                AllRangeFiles.append(rngfile)
    AllFilesWritten = list(AllFilesWritten)
    
    # Account for samkit run which produces no output files.  We still want
    # one info file written, using the name of the original input file to
    # generate an output file name.  In this case, every stream object will 
    # have an empty set for its "fileswritten" property so we simply create
    # a "fake" entry for the first stream in opts.orderedStreams
    if len(AllFilesWritten) == 0:
        filesWritten = [os.path.join(outputdir, opts.inputfile)]
        opts.orderedStreams[0].fileswritten = set(filesWritten)

    # Add info entries from other input files (.ini, .rngf, .rngpf etc)
    opts.infofile.configureExtraFilesSection(AllRangeFiles)

    # 4.14.4 - Add an info file counter
    infoCounter = 0
    infoFilesWritten = set()
    for stream in opts.orderedStreams:
        if not stream.op(OP_NOUT):
            if stream.op(OP_INFO) or not stream.op(OP_NOINFO):
                for file in tuple(stream.fileswritten) + tuple(extraFiles):
                    # MP runs fill in too many file names (both fastq files for example
                    # so we address that here
                    file = re.sub(r'\.pp\.[12]\.fastq', '.fastq', file)
                    file = re.sub(r'\.[12]\.fastq', '.fastq', file)
                    ifile  = '%s.info' % file
                    if file.startswith(outputdir) and not os.path.exists(ifile):
                        infoCounter += 1
                        infoFilesWritten.add(ifile)
    for ifile in list(infoFilesWritten):
        infofile = opts.infofile.copy()
        infofile.stats['runtype'] = 'complete'
        infofile.configureRunSection(opts,
                                     outputdir,
                                     AllFilesWritten,
                                     AllRangeFiles)
        infofile.configureHeader(outputdir,
                                 stream.OutputFiles, 
                                 stream.count)
        infofile.write(filename=ifile)

    # 4.14.4 - If any info files written, a single ".z.streams.csv" is required
    if infoCounter > 0:
        writeStreamCSV(opts, outputdir, streamSTATS)


def writeCSVfiles(opts, outputdir):
    """
    Write RNAME statistics to .CSV file(s), if requested.  Note that this 
    function returns a list of each CSV fie written
    """

    rnameStats    = compute_rname_stats(opts)
    csvExtensions = CSVExtensionList(opts.orderedStreams)

    filename = os.path.split(opts.inputfile)[1]
    fileroot,fileext = os.path.splitext(filename)
    if fileext == '.gz':
        fileroot,fileext = os.path.splitext(fileroot)

    filesWritten = []
    for csvExt in csvExtensions:
        rnameList     = UniqueRNAMEList(opts.orderedStreams, csvExt)
        csvfile = os.path.join(outputdir, opts.inputfile + '.' + csvExt)
        if os.path.exists(csvfile) and os.path.isfile(csvfile):
            continue
        if csvfile not in filesWritten:
            filesWritten.append(csvfile)
        fp = csv.writer(open(csvfile, 'w'), delimiter=',')

        suffixList    = SuffixList(opts.orderedStreams, csvExt)
        headers = ['RNAME', 'p.min', 'p.max']
        for suffix in suffixList:
            headers.append(re.sub('^\.', '', suffix))
        fp.writerow(headers)
 
        for rname in rnameList:
            # FLatten tuple/list (for RNAME pairs)
            row = []
            for item in rname:
                if type(item) in (list,tuple):
                    row.extend(item)
                else:
                    row.append(item)
            for suffix in suffixList:
                # RNAME pair
                if len(rname) == 2:
                    row.append(rnameStats[RnamePairAsKey(rname)][suffix])
                # Single RNAME
                else:
                    row.append(rnameStats[RnameAsKey(rname)][suffix])
            fp.writerow(row)
    return filesWritten


def writeStreamCSV(opts, outputdir, stats):
    """
    Write the ".z.streams.csv" file, a tab-separated list of suffixes, stream
    names and record counts.

    Introduced with version 4.14.4
    """
    myExtension = ".z.streams.csv"

    # Info file first @PR line will contain the correct timestamp
    pr = opts.infofile.PRlines[0].split(" ")[1:3]
    tstamp = "%s_%s" % (pr[0], pr[1].replace(":","."))
    
    # Assmemble a tuple of lines for the header (each one stored as a list)
    headerLines = ( ["# %s" % opts.inputfile,],
                    ["# %s" % outputdir.rstrip("/"),],
                    ["# %s" % tstamp,],
                  )

    # Create output file and write header lines followed by stats lines
    fileroot,fileext = os.path.splitext(opts.inputfile)
    csvfile = os.path.join(outputdir, fileroot + myExtension)
    fp = csv.writer(open(csvfile, 'w'), delimiter='\t')
    for line in headerLines:
        fp.writerow(line)
    for statline in stats:
       fp.writerow(statline)


def preProcessStreams(opts, heartbeat):

    if opts.pct is None or opts.pctLim is None:
        return

    # Set up a UUE stream for internal use
    optionString = 'test,UUE,fr,iL,0,%s,sam' % opts.pctLim['iLcalcmax']
    internal = SAMStream('internal', optionString)
    tokens = Tokenize(optionString)
    internal.parseOperatorTokens(tokens)

    total = 0
    limit = int(opts.pctLim['sample_size'])
    numlen = len("%d" % limit)
    strlen = 2*numlen+3
    print "Pre-processing to determine percentiles... %0*d/%0*d" % (numlen,
                                                                    0, 
                                                                    numlen,
                                                                    limit),
    
    for pair in ReadSAMFile(opts.inputfile):
        pair = internal.examinenext(pair, opts)
        if pair is None:
            total += 1
        print "%s" % (strlen*'\b'),
        print "%0*d/%0*d" % (numlen, total, numlen, limit),
        sys.stdout.flush()
        # Break out of processing loop if pre-process limit (default 100K) is
        # reached
        if total == int(opts.pctLim['sample_size']):
            break
    print "...done"

    # Compute and report percentile values
    SetPercentileValues(internal, opts, total)
    items = SortByValue(opts.pct)
    for k,v in items:
        print "   %s = %d" % (k, v)

    # Apply percentile values to streams
    try:
        for stream in opts.orderedStreams:
            stream.ops['pctile']  = opts.pct.copy()
            stream.ops['iLrange'] = stream.evaluatediLrange
    except:
        ErrorHandler()


def ProcessStreams(opts, order, callername, extraArgs=[], inifiles=[]):
    """
    Main driver function.  Accepts final opts which contains the required
    Stream class objects, as well as the input/output files etc.
    """
    # Validate input options, and update if necessary.  The opts object
    # will be modified in-place
    validateOptions(opts)

    # Generate a final dummy stream - a "catch-all" whose sole function is
    # to track any records not captured by any other stream.  This stream
    # produces no output and is not recorded in stats, other than to show
    # the output counts
    dummyArgs = 'unstreamed,AA,sam,nout'
    dummytokens = Tokenize(dummyArgs)
    dummyStream = SAMStream('n/a', dummyArgs)
    dummyStream.parseOperatorTokens(dummytokens)
    opts.orderedStreams.append(dummyStream)

    # Preprocess streams (or a subset of them).  This function may return
    # immediately if certain operators are not set
    heartbeat = RecordCounter(opts.inputfile)
    heartbeat.header()
    preProcessStreams(opts, heartbeat)

    # Set up controller and writer processes if multi-processing enabled
    if opts.nproc == 1:
        heartbeat.message("Processing on a single core", True)
    else:
        if opts.QLOW > (opts.QHIGH * 0.7):
            opts.QLOW = int(round(opts.QHIGH * 0.7))

        Controller = MPController(heartbeat, numProc=opts.nproc,
                                             chunkSize=opts.chunksize)
        Controller.start()
        # Use a simple or ordered writer, depending on option
        if opts.mp_unsorted:
            Writer = Concatenator(Controller.results, Controller.writerConn[1],
                                  opts.nproc, Controller.Counter, Controller.getpids())
        else:
            Writer = OrderedWriter(Controller.results, Controller.writerConn[1],
                                   opts.nproc, Controller.Counter)
        Writer.defibrillate(heartbeat)
        Writer.start()
        pids = Controller.getpids()
        pids.append(Writer.pid)
        write_kill_script(pids)

        #print "PROCESS ID Controller (self): ", os.getpid()
        #print "PROCESS ID processors:", Controller.getpids()
        #print "PROCESS ID Writer:", Writer.pid
        #print "NUM PROCESSES:", opts.nproc
        #print "QUEUE HIGH WATERMARK:", opts.QHIGH
        #print "QUEUE LOW WATERMARK:", opts.QLOW

    # Check and set limit value
    if opts.limit is not None:
        full_limit = opts.limit
    else:
        full_limit = 0

    # Set up InfoFile run parameters
    opts.infofile.openELN(sys.argv[1:] + extraArgs)
    for file in inifiles:
        opts.infofile.extraFiles.append(FileAndIdentifier(file))
    infofilecopy = opts.infofile.copy()

    # Set up output directories.  Output directory names will be stored as
    # an attribute 'outputDirectories' in opts
    makeDirectories(opts, callername)

    # Set up initial counter and update heartbeat
    total = 0
    heartbeat.update(force=True)

    # In Multi-processor mode, begin iterating through the input file
    # and converting pairs of lines into StreamProcessor callable task
    # objects before adding to the process queue
    # Iterate through records in input file.
    if opts.nproc > 1:
        Controller.Send( { 'orderedStreams' : opts.orderedStreams,
                         })
        paircount = 0
        chunkID = 0
        dataBloc = []
        pairBloc = []
        # Attempt open as a gzipped file, otherwise fall back
        try:
            fp = gzip.open(opts.inputfile, 'rb')
            fp._read_gzip_header()	# Throws an exception if not gzipped
            fp.rewind()			# Rewind on success
            fpIterator = fp.readlines	# Store file line iterator
        except IOError:
            fp.close()
            fp = open(opts.inputfile, 'r')
            fpIterator = fp.xreadlines	# Store file line iterator

        for line in fpIterator():
            if line.startswith('@'):
                continue
            pairBloc.append(line)
            paircount += 1
            if paircount % 2 == 0:
                dataBloc.append(pairBloc)
                pairBloc = []
                total += 1
            if len(dataBloc) == opts.chunksize:
                Controller.add(StreamProcessor(chunkID, dataBloc, opts),
                               opts.QHIGH, opts.QLOW)
                dataBloc = []
                chunkID += 1
            # Break out of processing loop if --limit records processed
            if full_limit > 0 and full_limit == total:
                break
        fp.close()

        # Any remaining pairs go into a final task construct
        if len(dataBloc) > 0:
            Controller.add(StreamProcessor(chunkID, dataBloc, opts),
                           opts.QHIGH, opts.QLOW)
            dataBloc = []
            chunkID += 1
            total += len(dataBloc)

        Controller.finishQueue()
        Controller.wait()
        dataDict = Controller.Recv()
        Controller.finishProcesses()
        opts.orderedStreams = dataDict['orderedStreams']

    # In single-processor mode use HTSeq to iterate through the input file
    # returning record pairs for processing.  A tally is kept of record pairs
    # which do not get filtered out by any stream
    else:
        for pair in ReadSAMFile(opts.inputfile):
            total += 1
            for stream in opts.orderedStreams:
                pair = stream.next(pair, opts)
                if pair is None:
                    break
            heartbeat += 1
            heartbeat.update()

            # Break out of processing loop if --limit records processed
            if full_limit == heartbeat.count:
                break

    # Now write the full run output
    WriteOutputFiles(opts, heartbeat, total, finalOutput=True)

    # Clean up empty output directories
    cleanDirectories(opts)

    if opts.nproc > 1:
        remove_kill_script()


def WriteOutputFiles(opts, heartbeat, total, finalOutput=False):
    """Update heartbeat information, and write output files"""

    # Update heartbeat display...
    if finalOutput:
        heartbeat.end()
    else:
        heartbeat.update(force=True)

    # Make backup copy infofile object
    infofile = opts.infofile.copy()

    # Compute statistics and write .info files and .csvN files
    outputdir = opts.outputDirectories.pop()
    streamSTATS = reportStats(opts, total, heartbeat)
    csvfiles = writeCSVfiles(opts, outputdir)
    writeInfoFiles(opts, outputdir, finalOutput, streamSTATS, csvfiles)

    print "FULL RUN COMPLETE:",   outputdir

