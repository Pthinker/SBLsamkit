####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####


"""
Classes and functions specifically for multi-processing
"""

import os
import sys
import copy
import gzip
import time
import types
import random
import datetime
import StringIO
from Queue import Full as qFull
from Operators import *
from SAM import makeSAMpairFromStringTuple
from Stream import SAMStream
from Utils import AvailableRAM, EmptyList
from multiprocessing import Process, JoinableQueue, Queue, Pipe, Value, Lock, cpu_count

import cProfile

# Define safety limits for Queue sizes to prevent overflows or running
# out of RAM.  A very loose ballpark estimate stands at perhaps 1.4MB per
# thousand tasks waiting in the queue. 
#
# A system with ~50GB of ram supports a queue size of 500000 with ease.  Systems
# with less RAM aren't so well behaved.  We'll set the watermark to 10000 (10k)
# times the number of available RAM in GB (rounded).  The low water mark will
# be 50% of this value
WATERMARK_BASE = 10000
QHIGH_WATERMARK = 500000
QLOW_WATERMARK  = 350000

class MPController(object):
    """
    Main MP object which maintains various queue and task objects and info.
    Launches generic client task processors as independant processes, populates
    the input queue and manages the writer process.
    """
    def __init__(self, heartbeat=None, numProc=None, chunkSize=1):
        self.heartbeat = heartbeat	# Heartbeat output manager
        self.tasks = JoinableQueue()	# Tasks for processing to be added here
        self.results = Queue()		# Processed results accumulate here
        self.writerConn = Pipe()	# Direct pipe to writer manager
        self.writer = self.writerConn[0]	
        self.resultDict = {}		# Info and orderedStream objects
        if numProc is None:
            numProc = max(cpu_count() - 2, 1)
        self.nproc = numProc
        self.chunkSize = chunkSize
        if self.heartbeat:
            #heartbeat.Lock = Lock()
            heartbeat.message("Launching %d sub-processes" % self.nproc, True)
            heartbeat.message("Tasks will process pairs in chunks of %d" % self.chunkSize, True)
        else:
            print "CONTROLLER WILL LAUNCH %d sub-processes" % self.nproc
            print "Tasks will process pairs in chunks of %d" % self.chunkSize
            sys.stdout.flush()
        #self.configureQueueLimits()
        self.tasksRunning = 0
        self.Counter = Value('L', 0)	# Shared value: counts processed pairs
        self.recordsProcessed = 0
        tLock = Lock()
        tValue = Value('L', 0)
        self.workers = [ TaskProcessor(self.tasks, self.results, tValue, tLock) 
                                               for i in xrange(self.nproc) ]
             

    def configureQueueLimits(self):
        availGB = AvailableRAM()
        self.HIGH_WATERMARK = int(WATERMARK_BASE * availGB)
        self.LOW_WATERMARK = int(round(self.HIGH_WATERMARK * 0.50))

    def add(self, task, qhigh, qlow):
        """Block if Qsize above a certain high watermark in item size, and
           don't release until it has fallen below the low watermark"""
        try:
            qlen = self.tasks.qsize()
            if qlen > qhigh:
                print "Throttling input, reached HWM:", qhigh
                while qlen > qlow:
                    delay = random.randint(1,10)
                    time.sleep(delay)
                    qlen = self.tasks.qsize()
                print "Throttling released, down to LWM:", qlow
        except NotImplementedError:
            # Skip on Mac OS X (WARNING - use on OS X in testing only, queue 
            # size will max out at a paltry 32768 items)
            pass
        try:
            self.tasks.put(task)
            self.recordsProcessed += task.datalen
        except qFull:
            # While testing: we shouldn't hopefully end up here...
            print "ERR: queue full"
            sys.exit(-1)

    def finishQueue(self):
        for i in xrange(self.nproc):
            self.tasks.put(None)
        if self.heartbeat is not None:
            self.heartbeat.total = self.recordsProcessed
            newCount = self.recordsProcessed
            self.heartbeat.message("Definitive pair count: %d" % newCount, True)

    def start(self):
        i = 0
        for worker in self.workers:
            worker.start()
            i += 1
        self.tasksRunning = self.nproc

    def OneTaskDone(self):
        self.tasksRunning -= 1

    def wait(self):
        self.tasks.join()
        if self.heartbeat:
            self.heartbeat.count = self.recordsProcessed

    def finishProcesses(self):
        self.writer.close()

    def getpids(self):
        pids = []
        for w in self.workers:
            pids.append(w.pid)
        return pids

    def Send(self, obj):
        self.writer.send(obj)

    def Recv(self):
        return self.writer.recv()


class TaskProcessor(Process):
    """Generic task processing object. A number of these are launched by the
    MPController as task processes.  Each one grabs task objects from the
    input Queue, executes them, and places the results in the result queue
    for the writer to take.
    """

    def __init__(self, task_queue, results_queue, orderID, taskLock):
        Process.__init__(self)
        self.taskQ = task_queue
        self.resultQ = results_queue
        self.orderID = orderID
        self.taskLock = taskLock

    def run(self):
        while True:
            task = self.taskQ.get()
            if task is None:
                self.resultQ.put(None)
                self.taskQ.task_done()
                break
            result = task(self.pid)	# Tasks = callable StreamProcessor 
					# objects: objects which can be called
					# like a function.
            self.addResult(result)
            self.taskQ.task_done()
        return

    def addResult(self, result):
        resultID = result[0]
        self.resultQ.put(result)


class StreamProcessor(object):
    """This is the processing guts object.  It accepts a tuple of strings
    representing a record pair as read from the input file (by the MPController)
    and converts it in a SAMPair object.  It uses the orderedStreams to
    process the SAMPair just as the SP process would, except that results
    are stored in a dict and include the task ID, stream name, output data,
    stats amd other accumultaed information.  The output data will be pre-
    formatted so the OrderedWriter does not have to manipulate them.
    """

    def __init__(self, ID, data, options):
        """Initialize object with task ID (also known as the chunkID) and
        the data (a string tuple), along with the defined options"""
        self.ID = ID
        self.data = data
        self.datalen = len(data)
        self.options = options

        # If the --unsorted option is used, each process writes out its own
        # set of files, tagss with its idividual PID.  
        self.writeToFiles = options.mp_unsorted

    def __call__(self, PID):
        """The meat of the parallel-processing approach, a __call__ method
        which makes the task object callable.  It converts the input tuple
        into a SAMPair object, runs it through each stream for processing
        in the normal fashion, and then passes the results to the output queue.
        """
        i = 0
        pairs = 0
        outputdata = []
        for recordpair in self.data:
            pair = makeSAMpairFromStringTuple(recordpair, reorder=False)
            for stream in self.options.orderedStreams:
                # In SP mode, stream.next() returns a pair or None.  In MP
                # it's more complicated, we pass back an array of dicts where
                # each one deinfes a pair (or not) depending on whether it is 
                # filtered out by the stream.
                result = stream.next(pair, self.options)
                if result['matched']:
                    if stream.op(OP_NOUT):
                        continue

                    # Copy stats for passing back.
                    copy_of_stats = copy.deepcopy(stream.stats)
                    copy_of_global = copy.deepcopy(self.options.orderedStreams[0].globalstats)

                    # Reset original stats.  Each subset of stats will
                    # be integrated separately
                    EmptyList(stream.stats)
                    EmptyList(self.options.orderedStreams[0].globalstats)

                    # First handle FASTQ output
                    dataBucketFASTQ = []

                    # Store root filename
                    froot = result['output'][0]

                    if stream.op(OP_FASTQ) or stream.op(OP_FASTQPP):
                        if stream.op(OP_FASTQ):
                            newpair,froot = self.ProcessPair(OP_FASTQ, stream, froot, pair)
                        else:
                            newpair,froot = self.ProcessPair(OP_FASTQPP, stream, froot, pair)
                        if self.writeToFiles:
                            if stream.op(OP_FASTQ) and stream.op(OP_SH):
                                outputf1 = "%s.sh.fastq.PID.%d" %(froot,PID)
                                if not stream.op(OP_INFO):
                                    dataBucketFASTQ = [open(outputf1, "a"),
                                                       None,
                                                      ]
                                else:
                                    dataBucketFASTQ = [None,
                                                       None,
                                                      ]
                            elif stream.op(OP_FASTQPP):
                                outputf1 = "%s.pp.1.fastq.PID.%d" %(froot,PID)
                                outputf2 = "%s.pp.2.fastq.PID.%d" %(froot,PID)
                                if not stream.op(OP_INFO):
                                    dataBucketFASTQ = [open(outputf1, "a"),
                                                       open(outputf2, "a"),
                                                      ]
                                else:
                                    dataBucketFASTQ = [None,
                                                       None,
                                                      ]
                            elif stream.op(OP_FASTQ):
                                outputf1 = "%s.1.fastq.PID.%d" %(froot,PID)
                                outputf2 = "%s.2.fastq.PID.%d" %(froot,PID)
                                if not stream.op(OP_INFO):
                                    dataBucketFASTQ = [open(outputf1, "a"),
                                                       open(outputf2, "a"),
                                                      ]
                                else:
                                    dataBucketFASTQ = [None,
                                                       None,
                                                      ]
                        else:
                            if not stream.op(OP_INFO):
                                dataBucketFASTQ = [StringIO.StringIO(), 
                                                   StringIO.StringIO(),
                                                  ]
                            else:
                                dataBucketFASTQ = [None,
                                                   None,
                                                  ]
                        if not stream.op(OP_INFO):
                            newpair.writeFASTQ(dataBucketFASTQ, closeWhenDone=False)


                    # Now Handle SAM output
                    dataBucketSAM = []

                    if stream.op(OP_SAM) or stream.op(OP_SAMPP):
                        if stream.op(OP_SAM):
                            newpair,froot = self.ProcessPair(OP_SAM, stream, froot, pair)
                        else:
                            newpair,froot = self.ProcessPair(OP_SAMPP, stream, froot, pair)
                        if self.writeToFiles:
                            if stream.op(OP_SAMPP):
                                outputf = "%s.pp.sam.PID.%d" %(froot,PID)
                                if not stream.op(OP_INFO):
                                    dataBucketSAM = [open(outputf, "a"),]
                                else:
                                    dataBucketSAM = [None,]
                            # OP_SAM (no OP_PP)
                            else:
                                outputf = "%s.sam.PID.%d" %(froot,PID)
                                if not stream.op(OP_INFO):
                                    dataBucketSAM = [open(outputf, "a"),]
                                else:
                                    dataBucketSAM = [None,]
                        else:
                            if not stream.op(OP_INFO):
                                dataBucketSAM = [StringIO.StringIO(),]
                            else:
                                dataBucketSAM = [None,]
                        if not stream.op(OP_INFO):
                            newpair.writeSAM(dataBucketSAM[0], closeWhenDone=False)


                    result['output'][0] = froot
                    # Return results
                    if stream.op(OP_SAM) or stream.op(OP_SAMPP) or \
                       stream.op(OP_FASTQ) or stream.op(OP_FASTQPP):
                        if self.writeToFiles:
                            if stream.op(OP_INFO):
                                files_for_output = []
                            else:
                                files_for_output = result['output']
                            outputdata.append({ 'datastrings' : '',
                                                'files': files_for_output,
                                                'name': result['name'],
                                                'stats': copy_of_stats,
                                                'gzipped' : stream.op(OP_GZ),
                                                'sam,pp' : stream.op(OP_SAMPP),
                                                'fastq,pp' : stream.op(OP_FASTQPP),
                                                'sh' : stream.op(OP_SH),
                                                'globalstats': copy_of_global,
                                              })
                        else:
                            pairvalueList = []
                            for db in dataBucketFASTQ + dataBucketSAM:
                                if db is None:
                                    pairvalueList.append(None)
                                else:
                                    # If a StringIO object has nothing written 
                                    # to it, the getvalue() call will throw an 
                                    # exception about the object not having a 
                                    # buf attribute. In this case we append None
                                    try:
                                        vv = db.getvalue()
                                        pairvalueList.append(vv)
                                    except:
                                        pairvalueList.append(None)

                            # "info" operator quashes SAM,FASTQ output
                            if stream.op(OP_INFO):
                                pairvalueList = []
                                files_for_output = []
                            else:
                                files_for_output = result['output']
                            outputdata.append({ 'datastrings' : pairvalueList,
                                                'files': files_for_output,
                                                'name': result['name'],
                                                'stats': copy_of_stats,
                                                'gzipped' : stream.op(OP_GZ),
                                                'sam,pp' : stream.op(OP_SAMPP),
                                                'fastq,pp' : stream.op(OP_FASTQPP),
                                                'sh' : stream.op(OP_SH),
                                                'globalstats': copy_of_global,
                                              })

                    for db in dataBucketFASTQ + dataBucketSAM:
                        try:
                            db.close()
                        except:
                            pass

                    if not stream.op(OP_PASS):
                        break
        

        # No matching data.  We'll return an "empty" output dict
        if len(outputdata) == 0:
            stream = self.options.orderedStreams[0]
            empty = SAMStream('none', '')
            outputdata = [{ 'datastrings' : '',
                            'files': [],
                            'name': empty.name,
                            'stats': empty.stats,
                            'gzipped' : False,
                            'sam,pp' : False,
                            'fastq,pp' : False,
                            'sh' : False,
                            'globalstats': stream.globalstats
                      },]
        return self.ID, outputdata

    def ProcessPair(self, outputoperator, stream, fileroot, pair):
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
                    paircopy = paircopy.proper_pair(OP_SAM, stream.ops)
                except ValueError:
                    fileroot += '.__stampy_PP__'
                except:
                    raise
            elif outputoperator == OP_FASTQPP:
                paircopy = paircopy.proper_pair(OP_FASTQ, stream.ops)
            elif outputoperator == OP_FASTQ:
                paircopy = paircopy.reverse_complement()

        # Convert to ordered pair.
        for eachOp in self.orderOp(stream):
            format, op = eachOp
            if stream.op(op) and format == outputoperator:
                fileroot = '%s.ord.%s' % (fileroot, op.split(':')[1])
                paircopy = paircopy.ordered_pair(stream, op)
                break

        return paircopy, fileroot

    def orderOp(self, stream):
        """
        Returns a list of tuples containing each order operator (if any) and
        the sub-stream which it will operator on
        """
        ops = []
        if stream.op(OP_ORDER):
            for op in stream.ops['ops']:
                if op.startswith(OP_ORDER):
                    ops.append( (self.leftStream(stream, op), op) )
        return ops

    def leftStream(self, stream, operator):
        """
        Return the substream declared to the left of a specific operator.  So
        for instance if an operator sequence looks like "sam,order,LL,fastq"
        and we passed "order:LL" as the operator this should return "sam".  It's
        a way of letting us apply an operator to a single output sub-stream.
        Note that if no sub-stream was declared to the left of the provided
        operator it's assumed to be "sam"
        """
        substreams = (OP_SAM, OP_SAMPP, OP_FASTQ, OP_FASTQPP)
        idx = stream.ops['ops'].index(operator)
        for i in range(idx,0,-1):
            if stream.ops['ops'][i] in substreams:
                return stream.ops['ops'][i]
        return OP_SAM


class BasicWriter(Process):
    """
    Base writer class. Handles some common interations (updating stats,
    sending stats back to controller, etc.

    This is an INCOMPLETE parent class, designed to be subclassed.
    """
    def __init__(self, output_queue, controllerPipe, nproc, Counter):
        Process.__init__(self)
        self.resultQ = output_queue
        self.controller = controllerPipe
        self.clients = nproc
        self.Counter = Counter
        self.outputfiles = set()
        self.filehandles = {'sam': {}, 'fastq': {}}
        self.heartbeat = None

    def defibrillate(self, heartbeat):
        """
        Start the heartbeat.  Geddit?  Defibrillate?
        """
        self.heartbeat = heartbeat

    def updateObjectsToController(self):
        """Send updated objects (in a dict) back to the controller"""
        newobj = { 'orderedStreams' : self.orderedStreams,
                   'outputfiles' : list(self.outputfiles),
                 }

        self.controller.send(newobj)
        self.controller.close()

    def integrateStats(self, data):
        """Incorporate the stats from a stream object with basically a single
        completed pair analyzed into the complete set of Stream objects that
        we can ten transmit back to the controller process.
        """

        typeTags = ('AA', 'UU', 'UUE', 'UUD', 'UX', 'UnM', 'UrM',
                    'nMnM', 'rMrM', 'nMrM', 'nMA', 'rMA', 'XX',
                    'fr', 'rf', 'ff', 'rr')

        for pDict in data:
            for i in xrange(len(self.orderedStreams)):
                if self.orderedStreams[i].name == pDict['name']:
                    stream = self.orderedStreams[i]
                    # insert size stats
                    for sizePair in pDict['stats']['sizes']:
                        stream.stats['sizes'].append(sizePair)
                    # Rname stats
                    for statsKey in ('rnames', 'RNsingle', 'RNpairs'):
                        for key,value in pDict['stats'][statsKey].items():
                            try:
                                stream.stats[statsKey][key] += value
                            except KeyError:
                                stream.stats[statsKey][key] = value
                    # LL and ss match stats
                    for tag in typeTags:
                        try:
                            stream.stats[tag] += pDict['stats'][tag]
                        except KeyError:
                            stream.stats[tag] = pDict['stats'][tag]
                    # General counts and output files (var: fns)
                    stream.count += 1
                    fns = []
                    for filename in pDict['files']:
                        if stream.op(OP_SAM):
                            fns.append('%s.%s' % (filename, OP_SAM))
                        elif stream.op(OP_SAMPP):
                            fns.append('%s.pp.%s' % (filename, OP_SAM))
                        if stream.op(OP_FASTQPP):
                            if stream.op(OP_SH):
                                fns.append('%s.pp.sh.%s' % (filename, OP_FASTQ))
                            else:
                                fns.append('%s.pp.1.%s' % (filename, OP_FASTQ))
                                fns.append('%s.pp.2.%s' % (filename, OP_FASTQ))
                        elif stream.op(OP_FASTQ):
                            if stream.op(OP_SH):
                                fns.append('%s.sh.%s' % (filename, OP_FASTQ))
                            else:
                                fns.append('%s.1.%s' % (filename, OP_FASTQ))
                                fns.append('%s.2.%s' % (filename, OP_FASTQ))
                        for fname in fns:
                            stream.fileswritten.add(fname)
                            stream.outputfilenames.add(fname)
                    break

        # Now add the so-called global stats from each dict - we only need
        # update the globalstats list in the first stream object
        for pDict in data:
            stream = self.orderedStreams[0]
            for tag in typeTags:
                try:
                    stream.globalstats[tag] += pDict['globalstats'][tag]
                except KeyError:
                    stream.globalstats[tag] = pDict['globalstats'][tag]

        # Update the heartbeat object with the number of record pairs
        # written. This will trigger a console update automatically when
        # needed.
        self.heartbeat.count = stream.globalstats['AA']
        self.heartbeat.update()


class Concatenator(BasicWriter):
    """
    Faux-writer class.  It processes record data (which it should be noted does
    NOT contain any actual pair data) for stats and updates stats in the parent
    controller process.  Also concatenates all files from each registered 
    processor on close.
    """
    def __init__(self, output_queue, controllerPipe, nproc, Counter, pidList):
        BasicWriter.__init__(self, output_queue, controllerPipe, nproc, Counter)
        self.pidList = pidList

    def run(self):
        """Begin processing results from the output queue.  Blocks if the
        queue is empty but client processes are still running.  Client counter
        is decremented any time a "poison pill" (None object) is retrieved from
        the results queue, until the client counter goes to zero.
        Each result comes in the form of a tuple containing the task ID and
        data as a list of dicts.  The task ID is ingored (used by the
        OrderedWriter class). No data is written, stats are simply updated.
        (which may not match the original input order).
        """
        # Get data objects (in a dict) from the controller process 
        dataDict = self.controller.recv()
        self.orderedStreams = dataDict['orderedStreams']

        ID = None
        data = None
        output_compressed = set()
        output_normal     = set()
        while self.clients:
            result = self.resultQ.get()
            if result is None:
                self.clients -= 1
                continue
            ID, data = result
            for pDict in data:
                if pDict['gzipped']:
                    for filename in pDict['files']:
                        output_compressed.add(filename)
                else:
                    for filename in pDict['files']:
                        output_normal.add(filename)
                for filename in pDict['files']:
                    self.outputfiles.add(filename)
            
            self.integrateStats(data)
            self.Counter.value += len(data)

        # Now concatenate any output files together
        if self.heartbeat is not None:
            self.heartbeat.message("Beginning file block merging..", True)

        fcount = 0
        blkavg = 0
        for extension in ('sam', 
                          'pp.sam',
                          '1.fastq',
                          '2.fastq',
                          'pp.1.fastq',
                          'pp.2.fastq',
                          'sh.fastq',
                          'sh.pp.fastq'):
            fc,ba = self.concatenate(output_compressed, extension, do_gzip=True)
            fcount += fc
            blkavg += ba
            fc,ba = self.concatenate(output_normal, extension, do_gzip=False)
            fcount += fc
            blkavg += ba

        if self.heartbeat is not None and fcount > 0:
            self.heartbeat.message(
                    "Merged %d blocks (avg) in each of %d output files" % 
                    (int(round(blkavg * 1.0 / fcount)), fcount), True)
        

        # Send updated data (stats mainly) via the pipe directly back to
        # the MPController object, close filehandles and finish up.
        self.updateObjectsToController()

    def concatenate(self, output_files, extension='sam', do_gzip=False):
        fcount = 0
        blkavg = 0
        for filename in output_files:
            ofile = '%s.%s' % (filename, extension)
            try:
                if os.path.isfile(ofile):
                    os.unlink(ofile)
            except OSError:
                print "Unable to write '%s', skipping..." % ofile
                continue
            fList = []
            command = ""
            for PID in self.pidList:
                subfile = "%s.PID.%d" % (ofile,  PID)
                if os.path.isfile(subfile):
                    fList.append(subfile)
            if len(fList) > 0:
                fcount += 1
                blkavg += len(fList)
                if do_gzip:
                    os.system("cat %s | gzip -9 > %s.gz" % (' '.join(fList),
                                                            ofile))
                else:
                    print "cat %s > %s" % (' '.join(fList), ofile)
                    os.system("cat %s > %s" % (' '.join(fList), ofile))
                for f in fList:
                    try:
                        os.unlink(f)
                    except:
                        pass
        return (fcount, blkavg)
                

class SimpleWriter(BasicWriter):
    """
    Writer class for the output process.  It takes record pairs from the
    output queue (the StreamProcessor tasks place them there in a dict which
    includes a unique ID and a status result).  The output data may be an
    empty array (no match, nothing to output).  Items are writtein in the order
    in which they arrive in the queue.
    """

    def run(self):
        """Begin processing results from the output queue.  Blocks if the
        queue is empty but client processes are still running.  Client counter
        is decremented any time a "poison pill" (None object) is retrieved from
        the results queue, until the client counter goes to zero.
        Each result comes in the form of a tuple containing the task ID and
        data as a list of dicts.  The task ID is ingored (used by the
        OrderedWriter class) and the data is writtein in the order it appears
        (which may not match the original input order).
        """
        # Get data objects (in a dict) from the controller process 
        dataDict = self.controller.recv()
        self.orderedStreams = dataDict['orderedStreams']

        ID = None
        data = None
        while self.clients:
            result = self.resultQ.get()
            if result is None:
                self.clients -= 1
                continue
            ID, data = result
            # Data sequence is unimportant, simply write it out and proceed
            self.writePairs(data)

        # Send updated data (stats mainly) via the pipe directly back to
        # the MPController object, close filehandles and finish up.
        self.updateObjectsToController()
        self.closeFileHandles()

    def writePairs(self, data):
        # Pairs can be 1 (SAM only), 2 (FASTQ only) or 3 (FASTQ and SAM). In
        # the thrd case FASTQ will be processed first)
        for pDict in data:
            for filename in pDict['files']:
                self.outputfiles.add(filename)
                self.writePairSAM(pDict, filename)
                self.writePairFASTQ(pDict, filename)
        self.integrateStats(data)
        self.Counter.value += len(data)

    def writePairSAM(self, pDict, fileroot):
        if len(pDict['datastrings']) == 0:
            return
        # See writePairs above - the SAM record is always the last one
        if pDict['sam,pp']:
            fileroot = '%s.pp' % fileroot
        if fileroot not in self.filehandles['sam'].keys():
            fn = '%s.%s' % (fileroot, 'sam')
            if pDict['gzipped']:
                fp = gzip.open('%s.gz' % fn, 'wb')
            else:
                fp = open(fn, 'w')
            self.filehandles['sam'][fileroot] = fp
        else:
            fp = self.filehandles['sam'][fileroot]
        if pDict['datastrings'][-1] is not None:
            print >> fp, pDict['datastrings'][-1],

    def writePairFASTQ(self, pDict, fileroot):
        if len(pDict['datastrings']) < 2:
            return
        # See writePairs above - the FASTQ records only exist if there
        # are at least 2 elements in pDict['datastrings']
        if pDict['sh']:
            fileroot = '%s.sh' % fileroot
        if pDict['fastq,pp']:
            fileroot = '%s.pp' % fileroot
        if fileroot not in self.filehandles['fastq'].keys():
            if pDict['sh']:
                fn1 = '%s.%s' % (fileroot, 'fastq')
                fn2 = ''
                if pDict['gzipped']:
                    fp1 = gzip.open('%s.gz' % fn1, 'wb')
                    fp2 = fp1
                else:
                    fp1 = open(fn1, 'w')
                    fp2 = fp1
            else:
                fn1 = '%s.%s' % (fileroot, '1.fastq')
                fn2 = '%s.%s' % (fileroot, '2.fastq')
                if pDict['gzipped']:
                    fp1 = gzip.open('%s.gz' % fn1, 'wb')
                    fp2 = gzip.open('%s.gz' % fn2, 'wb')
                else:
                    fp1 = open(fn1, 'w')
                    fp2 = open(fn2, 'w')
            self.filehandles['fastq'][fileroot] = [fp1, fp2]
        else:
            fp1, fp2 = self.filehandles['fastq'][fileroot]
        if pDict['datastrings'][0] is not None and \
           pDict['datastrings'][1] is not None:
            print >> fp1, pDict['datastrings'][0],
            print >> fp2, pDict['datastrings'][1],

    def closeFileHandles(self):
        for fp in self.filehandles['sam'].values():
            if fp is not None:
                fp.close()
        for fp in self.filehandles['fastq'].values():
            for handle in fp:
                # Both handles may be the same (if shuffled output via OP_SH is
                # enabled), so catch any exception on attempting to close the
                # samepl file twice
                try:
                    handle.close()
                except:
                    pass


class OrderedWriter(SimpleWriter):
    """
    Writer class for the output process.  It takes record pairs from the
    output queue (the StreamProcessor tasks place them there in a dict which
    includes a unique ID and a status result).  The output data may be an
    empty array (no match, nothing to output).  Items arriving out of order
    are identified by a missing value in the sequence of unique IDs and stored
    in an output cache.  This cache is then drained as far as possible each
    time a new in-sequence result is taken from the output queue.
    """

    def __init__(self, output_queue, controllerPipe, nproc, Counter):
        SimpleWriter.__init__(self, output_queue, controllerPipe, nproc,Counter)
        self.cache = []
        self.cacheIDX = []

    def run(self):
        """Begin processing results from the output queue.  Blocks if the
        queue is empty but client processes are still running.  Client counter
        is decremented any time a "poison pill" (None object) is retrieved from
        the results queue, until the client counter goes to zero.
        Each result comes in the form of a tuple containing the task ID and
        data as a list of dicts.  If the task ID is in sequence it is written
        immediately.  If out of sequence it is stored in a cache (basically
        a local list of objects).  Each out-of-sequence item is stored in the
        cache until the next in-sequence object arrives.  After it is written
        attempts are made to drain the cache as far as possible by locating
        each in-sequence ID.
        """

        # Get data objects (in a dict) from the controller process 
        dataDict = self.controller.recv()
        self.orderedStreams = dataDict['orderedStreams']
        orderedID = 0

# TO DO - perhaps treat empty results differently?  Place the ID in a "discard"
# list and do the stats update immediately.  Then when incrementing the 
# orderedID, if the updated value is in the discard list, increment again..
        # Begin procesing results queue
        ID = None
        data = None
        c = 0
        while self.clients:
            result = self.resultQ.get()
            if result is None:
                self.clients -= 1
                continue
            c += 1
            ID, data = result
            while self.clients and ID != orderedID:
                result = self.resultQ.get()
                if result is None:
                    self.clients -= 1
                    continue
                c += 1
                self.cacheIDX.append(ID)
                self.cache.append(data)
                ID, data = result

            # Data is next in sequence, write it out and proceed
            self.writePairs(data)
            orderedID += 1
            while orderedID in self.cacheIDX:
                idx = self.cacheIDX.index(orderedID)
                self.cacheIDX.pop(idx)
                data = self.cache.pop(idx)
                self.writePairs(data)
                orderedID += 1

        # Processing is completed but the cache may not be empty.  Drain it
        # now (it should contain any missing objects at this point)
        if len(self.cacheIDX):
            while orderedID in self.cacheIDX:
                idx = self.cacheIDX.index(orderedID)
                self.cacheIDX.pop(idx)
                data = self.cache.pop(idx)
                self.writePairs(data)
                orderedID += 1

        # Send updated data (stats mainly) via the pipe directly back to
        # the MPController object, close filehandles and finish up.
        self.updateObjectsToController()
        self.closeFileHandles()



def write_kill_script(pids):
    """This is a temporary measure to provide a "kill script".  If the 
    application crahses, it will leave a series of child processes running
    in which case the "kill_script.sh" script generated here can be used
    to cleanup.  Once executed, the script will erase itself to avoid a
    potentially unpleasant side-effect where if runat some arbitrary point
    in the future, unexpected processes might be terminated.

    Eventaully we'll need this to be replaced with an auto-cleanup feature
    which for example might be called on any un-caught exception."""
    fp = open('kill_script.sh', 'w')
    for pid in pids:
        print >> fp, "kill -TERM %d" % pid
    print >> fp, "kill -TERM %d" % os.getpid()
    print >> fp, "rm -f kill_script.sh"
    fp.close()

def remove_kill_script():
    """Assuming all went well and no exceptions were raised, the kill script
    won't be required and can be removed
    """
    try:
        os.unlink('kill_script.sh')
    except:
        pass


