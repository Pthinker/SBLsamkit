####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
Miscellaneous utility functions and classes
"""

import re
import csv
import sys
import types
import os.path
import datetime
import string
import locale
import tokenize
from random import uniform
from numpy import inf
from StringIO import StringIO
from Operators import OP_CSV
from Exceptions import *


def Ellipsis(input, length=50):
    """
    Truncates the input string if longer than the specified length (default 50)
    to length-3 and appends an ellipsis 
    """
    if len(input) > length:
        return input[:length-3] + '...'
    else:
        return input

def Integer(value):
    """
    Integer parser which accepts standard multiplier suffixes and applies them
    to the numerical values.  So for example 100K or 100k would set the value
    to an integer 1000000.  Accepted multipliers, ignoring case, are:
    K (times 1000)
    M (times 1000000)
    NOTE:  Accepts floating-point notation, rounds to nearest int.  Thus 
    an input value 0.5M is valid (500000)
    """
    try:
        value = value.upper()
        if value.endswith('K'):
            value = float(value.replace('K','')) * 1000
        elif value.endswith('M'):
            value = float(value.replace('M','')) * 1000000
        else:
            value = float(value)
        value = int(value)
    except ValueError:
        raise
    return value

def LocalizedInteger(value):
    """Returns a string version of the integer with seperators"""
    locale.setlocale(locale.LC_ALL, "")
    return locale.format('%d', value, True)


def iLrangeAsString(iLrange):
    """Returns a pretty-print string representing the iLrange which is passed
    in as a list of two values.  String length is not guaranteed to be below
    any limit.
    """
    iLmin, iLmax = iLrange
    if type(iLmin) is types.IntType:
        if iLmin > 0 and iLmin % 1000 == 0:
            iLmin /= 1000
            iLmin = '%dk' % iLmin
        else:
            iLmin = '%d' % iLmin

    if type(iLmax) is types.FloatType:
        if iLmax == inf:
            iLmax = 'inf'
        else:
            iLmax = int(iLmax)
    if type(iLmax) is types.IntType:
        if iLmax > 0 and iLmax % 1000 == 0:
            iLmax /= 1000
            iLmax = '%dk' % iLmax
        else:
            iLmax = '%d' % iLmax

    return '%s-%s' % (iLmin, iLmax)


def ExpandPath(filename):
    """
    Returns the equivalent file or directory name with tilde expansion 
    completed.  Thus ~/foo.sam becomes /home/username/foo.sam
    """
    try:
        return os.path.expanduser(filename)
    except AttributeError:
        return filename


def FileExists(filename):
    """
    Test whether a file exists
    """
    return os.path.exists(filename)


def FileReadable(filename):
    """
    Test whether a file can be opened and its contents read
    """
    try:
        f = open(filename, 'r')
        f.close()
        return True
    except:
        return False


def Error(message, exitval=-1, options=None):
    """
    Given an input message (string), display the message with a suitable
    label, and exit with the provided exit value.  If an OptionParser 
    object is passed in, call its print_help() method before displaying
    the specific message provided.
    """
    print >> sys.stderr, ""
    if options:
        options.print_help() 
        print >> sys.stderr, ""
    print >> sys.stderr, "ERROR: %s" % message
    print >> sys.stderr, ""
    sys.exit(exitval)


def ErrorHandler():
    """
    A custom error handler function which can be called without arguments
    when any exception is caught.  Extracts and displays custom exception
    information if applicable, otherwise just re-raises the exception
    """
    exception = sys.exc_info()
    if str(exception).find('SBLsamkit.Exceptions') >= 0:
        Error(exception[1])
    else:
        raise


def EmptyList(L):
    """
    Empty a list of simple object types - may be ints, floats, strings,
    lists, tuples or dicts:
    int -> 0
    float -> 0.0
    string -> ''
    list -> []
    tuples -> ()
    dict -> {}
    The list is emptied in-place
    """
    for k in L.keys():
        if type(L[k]) is types.DictType:
            L[k] = {}
        elif type(L[k]) is types.ListType:
            L[k] = []
        elif type(L[k]) is types.TupleType:
            L[k] = ()
        elif type(L[k]) is types.StringType:
            L[k] = ''
        elif type(L[k]) is types.FloatType:
            L[k] = 0.0
        elif type(L[k]) is types.IntType:
            L[k] = 0


def EstimateRecordCount_SAM(filename):
    """
    A function to estimate the total number of records in the input file
    *without* loading it all into memory.  We do this by reading the
    first 2000 lines (ostensibly 1000 records) and dividing the total file
    size (in bytes) by the number of bytes we've read in over times 500.  It's 
    not absolutly accurate of course, but will suffice for an estimate.
    """
    # Read up to 2000 lines from file, skipping headers, and count the total 
    # number of bytes read.  We also record the file pointer position after 
    # skipping the headers, and will remove this value (plus the length of the 
    # non-header line) from the total file size to improve accuracy
    filehandle = open(filename, 'r')
    line = filehandle.readline()
    while line.startswith('@'):
        line = filehandle.readline()
    lines_read = 1
    bytes_read = len(line)
    fp_position = filehandle.tell() - bytes_read
    for i in range(0,2000):
        line = filehandle.readline()
        if len(line) == 0:
            break
        bytes_read += len(line)
        lines_read += 1
    filehandle.close()

    # Return the estimated record count
    if bytes_read == 0 or lines_read == 0:
        return 0
    else:
        lines_read /= 2
        return FileSize(filename, fp_position) * lines_read / bytes_read


def FileSize(filename, offset=0):
    """
    Given a filename, return the file size in bytes.
    """
    return os.path.getsize(filename) - offset


def RnameAsKey(rnameTuple):
    """
    Given an RNAME match tuple in the form RNAME, p_min, p_max, returns a
    standard string representation suitable as a hash key.  RNAME is a string
    and p_min and p_max are integers.  Any time a numpy.inf is encountered
    it is converted to -1
    """
    for i in range(1,3):
        if rnameTuple[i] == inf:
            rnameTuple[i] = -1
    return "%s<>%d<>%d" % (rnameTuple[0], rnameTuple[1], rnameTuple[2])


def RnamePairAsKey(tupleofRnameTuples):
    """
    Given a tuple of 2 RNAME match tuples (each in the form RNAME, p_min, p_max,
    returns a standard string representation suitable as a hash key, basically
    consisting of the RnameAsKay() output for each in the pair concatenated
    together.
    """
    keys = []
    for tup in tupleofRnameTuples:
        keys.append(RnameAsKey(tup))
    return '<>'.join(keys)


def KeyAsRname(keyString):
    """
    The inverse of the RnameAsKey function, this accepts a string in the 
    form "RNAME<>p_min<>p_max" and returns a valid RNAME, p_min, p_max tuple
    where RNAME is a string, p_min and p_max are integers.  Any time a "-1"
    or "inf" string is detected for p_max, it is converted to numpy.inf
    """
    elements = keyString.split('<>')
    for i in range(1,3):
        if elements[i] == 'inf':
            elements[i] = inf
        else:
            elements[i] = int(elements[i])
    return (elements[0], elements[1], elements[2])


def SuffixList(streamList, operator=None):
    """
    Constructs a set of unique suffixes defined in the provided list of streams.
    If the optional operator is provided, only those streams which include
    that operator string in their operator list will be included. 

    If a string which ends with 'csv' is passed in, only those stream whose
    ops['csvExt'] attribute match will be included.

    Note that passing operator None to stream.op() always returns True.
    """
    suffixList = []
    for stream in streamList:
        if stream.op(operator):
            suffix = stream.ops['suffix']
            if suffix not in suffixList:
                suffixList.append(suffix)
        elif operator.endswith(OP_CSV) and operator == stream.ops['csvExt']:
            suffix = stream.ops['suffix']
            if suffix not in suffixList:
                suffixList.append(suffix)
    return suffixList


def CSVExtensionList(streamList):
    """
    Constructs a set of unique CSV output file extensions defined in the 
    provided list of streams
    """
    csvExtensions = set()
    for stream in streamList:
        if stream.op(OP_CSV):
            csvExtensions.add(stream.ops['csvExt'])
    return csvExtensions


def UniqueRNAMEList(streamList, csvExtension=None):
    """
    Constructs a list of unique RNAME,p_min,p_max patterns.  If an optional
    csvExtension string is provided, only streams whose ops['csvExt'] attribute
    match will be included when examining each stream.
    """
    rnameList = []
    for stream in streamList:
        if csvExtension is not None and csvExtension != stream.ops['csvExt']:
            continue
        if len(stream.ops['rnames']) > 0:
            if stream.ops['RNpairs']:
                for i in range(0, len(stream.ops['rnames']),2):
                    rnitem0 = stream.ops['rnames'][i]
                    rnitem1 = stream.ops['rnames'][i+1]
                    rnitem = (rnitem0, rnitem1)
                    if rnitem not in rnameList:
                        rnameList.append(rnitem)
            else:
                for rnitem in stream.ops['rnames']:
                    if type(rnitem) in (tuple,):
                        rnames = rnitem
                    else:
                        rnames = (rnitem,)
                    for rn in rnames:
                        if rn not in rnameList:
                            rnameList.append(rn)
    return rnameList


def Tokenize(tokenstring, specialTokens=(':',)):
    """
    Splits incoming token string at commas.  Tokens with special characters
    included in the specialTokens tuple (options, default is a single ':' in
    a one-element tuple) are subsequently split again (e.g. the ':' between 
    range pairs). 
    """
    tokens = tokenstring.split(',')
    for specialToken in specialTokens:
        for i in range(0, len(tokens)):
            if specialToken in tokens[i]:
                subtoken = tokens.pop(i)
                subtokens = subtoken.split(specialToken)
                for j in range(0,len(subtokens) - 1):
                    subtokens[j] += specialToken
                for j in range(len(subtokens)-1, -1, -1):
                    tokens.insert(i, subtokens[j])
    return tokens


def IterateInPairs(inputList):
    """
    Generator function which accepts a flat list of tokens (strings or 
    whatever) and returns a tuple of item pairs, taking the first two items in 
    the list for each pass.  Thus an input list [1,2,3,4,5,6] in a foreach loop
    returns (1,2), then (3,4), then (5,6)
    """

    # Make a copy of the list otherwise we'll be modifying the input list!
    tokenList = list(inputList)
    
    while True:
        try:
            yield (tokenList.pop(0), tokenList.pop(0))
        except IndexError:
            raise StopIteration

 
def ProcessPercentileOption(opts):
    """Process --percentile option, supports setting variables
    only letters (upper and lower case), numbers, and underscores allowed)
    """
    allowed_characters = string.letters + string.digits + '_'
    if opts.pct:
        percentileTokens = Tokenize(opts.pct)
        if len(percentileTokens) < 4 or len(percentileTokens) % 2:
            raise PercentileParserException(
                      "Invalid --percentile declaration '%s'" % opts.pct)
        try:
            sample_sizeToken = percentileTokens.pop()
            sample_size = Integer(sample_sizeToken)
        except ValueError:
            raise PercentileParserException(
                      "Invalid integer value '%s' for percentile sample size" %
                      sample_sizeToken)
        try:
            iLcalcmaxToken = percentileTokens.pop()
            iLcalcmax      = Integer(iLcalcmaxToken)
        except ValueError:
            raise PercentileParserException(
                      "Invalid integer value '%s' for percentile iL calc max" %
                      iLcalcmaxToken)
       
        opts.pct = {}
        opts.pctLim      = { 'iLcalcmax'  : iLcalcmax,
                             'sample_size': sample_size }

        for (percentile, variable) in IterateInPairs(percentileTokens):
            for ch in variable:
                if ch not in allowed_characters:
                    raise PercentileParserException(
                          "Invalid character '%s' in percentile var '%s'" % (
                              ch, variable))
            try:
                opts.pct[variable] = float(percentile)
            except:
                raise PercentileParserException(
                      "Invalid float value '%s' for percentile var '%s'" % 
                      (percentile, variable))
    return opts


def MakeExpression(rawstring):
    """
    Given an input string, converts to a valid math expression where tokens
    (characters) which are not basic math symbols and are not numerical 
    values are treated as variables and replaced with a valid expression
    to reference those variables
    """

    mathTokens = ('+', '-', '*', '/', '(', ')')
    expression = ['',]
    generator = tokenize.generate_tokens(StringIO(rawstring).readline)
    for toknum,tokval,_,_,_ in generator:
        if toknum in (tokenize.NAME, tokenize.NUMBER):
            if expression[-1] in mathTokens:
                expression.append(tokval)
            else:
                expression[-1] = expression[-1] + tokval
        elif toknum == tokenize.OP and tokval in mathTokens:
            expression.append(tokval)
        elif toknum == tokenize.ENDMARKER:
            break
        else:
            raise ValueError('Invalid token "%s"' % tokval)
    return expression


def SortByValue(dictionary):
    """
    As the name suggests, sort a dictionary *numerically* by value, returning 
    a list of key,value pairs, just as dict.items() would do
    """
    items = dictionary.items()
    items.sort(cmp=lambda x,y: cmp(x[1], y[1]))
    return items


def FoundIn(needle, haystack, caseSensitive=False):
    """
    Return True if the needle (a string) occurs in haystack (a list of strings)
    OR the needle with a leading '^' character removed occurs in haystack
    """
    targets = []
    if not caseSensitive:
        search = needle.upper()
        for item in haystack:
            targets.append(item.upper())
    else:
        search = needle
        targets = haystack
    return search in targets or search.replace('^','') in targets


def TokensFromCSV(filename):
    """
    Reads a CSV file and returns the contents as a list of items.  Whether the
    file contents is split across multiple lines or not, a single list is
    always returned.  Lines beginning with the comment character 3 are skipped.
    Also returns the first skipped line (the file identifier) as a string, if
    it exists, or a blank string otherwise
    """
    csvtokens = []
    filehandle = open(filename, 'rb')

    # Store the first line if it begins with a #
    firstline = filehandle.readline().rstrip()
    filehandle.seek(0)
    if not firstline.startswith('#'):
        firstline = ''

    # Now skip each line that begins with a comment
    location = 0
    if firstline != '':
        skipline = '%s' % firstline
        while skipline.startswith('#'):
            location = filehandle.tell()
            skipline = filehandle.readline().rstrip()
            continue
    
    # Determine the input file dialect
    filehandle.seek(location)
    sample = filehandle.read(1024)
    dialect = csv.Sniffer().sniff(sample)
    # We could use the Sniffer.has_header() class method to determine if the
    # input file begins with a header row, e.g.:
    #     headerExists = csv.Sniffer().has_header(sample)
    # However the samkit4 spec calls for input files which may have one or more
    # header rows at the beginning,(detected if beginning with a # symbol) so 
    # instead we "seek" to the first non-header line before passing control to      # to the csv reader

    # Read input file using the determined dialect, append each line as as
    # list of tokens, beginning at the first non-header line
    filehandle.seek(location)
    reader = csv.reader(filehandle, dialect)
    for row in reader:
        csvtokens.extend(row)
    return csvtokens, firstline


def FileAndIdentifier(filename):
    """
    Returns tuple of filename and first line (ID) if a comment, or en empty
    string.  Thus it will return either (filename, '# blah blah') or 
    (filename, '')
    """
    try:
        fH = open(filename, 'r')
        id = fH.readline().rstrip()
        if not id.startswith('#'):
            id = ''
        close(fH)
    except:
        id = ''
    return (filename, id)


def AvailableRAM():
    """
    Returns the total available RAM in the system (not free RAM).  Returns
    the value from in /proc/meminfo in nearest GB (Linux), or a hardcoded
    value otherwise (for testing e.g. on OS X)
    """
    try:
        fp = open('/proc/meminfo')
    except IOError:
        return 8.0

    memGB = 0.0
    for line in fp.xreadlines():
        line = line.strip()
        if line.startswith('MemTotal:'):
            fp.close()
            values = line.split()
            if len(values) == 3:
                memGB = int(values[1])
                if values[2].lower() == 'kb':
                    memGB *= 1024
                if values[2].lower() == 'mb':
                    memGB *= 1024 * 1024
                if values[2].lower() == 'gb':
                    memGB *= 1024 * 1024 * 1024
                return round(memGB / (1024.0 * 1024.0 * 1024.0))
    # Horrible, but a default will suffice for the moment.
    fp.close()
    return 1.0


class RecordCounter(object):
    """
    This basic class tracks pairs read, and compares against total (or 
    estimated total), to report a progress "heartbeat" to the console
    """

    INTERVAL = 50000

    def __init__(self, filename, limit=None):
        self.count = 0
        self.total = EstimateRecordCount_SAM(filename)
        if limit is not None:
            self.total = limit * 1000
            self.limit = limit
        else:
            self.limit = None
        self.startTime,self.startTime_str = self.now
        self.endTime = None
        self.previousValue = 0
        self.Lock = None

    def header(self):
        self.lock()
        print "%s: BEGIN - Estimated pair count: %d" % (self.startTime_str,
                                                        self.total)
        sys.stdout.flush()
        self.unlock()

    def footer(self):
        self.lock()
        print "%s: END - Actual pair count: %d" % (self.now[1], self.count)
        sys.stdout.flush()
        self.unlock()

    def begin(self):
        self.header()
        self.update(force=True)

    def end(self):
        self.total = self.count
        self.update(force=True)
        self.endTime = self.now[0]
        self.footer()

    def message(self, text, timestamped=False):
        self.lock()
        if timestamped:
            print "%s: INFO - %s" % (self.now[1], text)
        else:
            print '  ** %s **' % text
        sys.stdout.flush()
        self.unlock()

    def update(self, force=False):
        """Update heartbeat display if requested"""
        # update may not be called with self.count exactly being a multiple
        # of self.INTERVAL.  So we check the number of times self.INTERVAL
        # has been passed, and force an update when it changes...
        self.lock()
        divisors = self.count / self.INTERVAL
        if divisors > self.previousValue:
            force = True
            self.previousValue = divisors
        if force is True or self.count % self.INTERVAL == 0:
            t_elapsed = self.time_elapsed()
            try:
                percent = min(100.0, (self.count*100.0)/self.total)
            except ZeroDivisionError:
                percent = 0.0
            if self.count >= 1000:
                processed = "%dK" % (self.count / 1000)
            else:
                processed = "%d" % self.count
            print '%10s records (%6.2f%%) - elapsed time: %02dm %02ds' % (
                                                           processed,
                                                           percent,
                                                           t_elapsed['minutes'],
                                                           t_elapsed['seconds'],
                                                          )
            sys.stdout.flush()
        self.unlock()

    def time_elapsed(self):
        """Returns the elapsed time in a dict of minutes and seconds"""
        interval = {}
        elapsed = self.now[0] - self.startTime
        interval['minutes'] = elapsed.seconds / 60
        interval['seconds'] = elapsed.seconds % 60
        return interval

    def __add__(self, value):
        """Increments the self.count variable.  This allows the RecordCounter
        object to be simply incremented, e.g.:
        r = RecordCounter()
        r += 1
        """
        self.count += value
        return self

    def timestr(self, timestamp):
        """Returns a human-redable string of the provided timestamp"""
        return timestamp.strftime('%d/%m/%Y %H:%M:%S')
        
    def complete(self):
        return self.count >= self.total

    def lock(self):
        """Acquire multiprocessing lock, if self.Lock has been set as a Lock"""
        if self.Lock:
            self.Lock.acquire()

    def unlock(self):
        """Release multiprocessing lock, if self.Lock has been set as a Lock"""
        if self.Lock:
            self.Lock.release()

    @property
    def now(self):
        """Returns a native timestamp, and a formatted string version of it
        in a tuple
        """
        timestamp = datetime.datetime.today()
        return timestamp, self.timestr(timestamp)

    @property
    def TotalRunTimeInSeconds(self):
        """This property returns the total run time in seconds, if the 
        self.endTime attribute has been set (which is done when self.end() is
        called).  If this attribute has not been set, the current elapsed time
        is returned.  In both cases, only the 'seconds' attribute of the 
        timedelta object is returned (an integer)
        """
        if self.endTime is None:
            elapsedTime = self.now[0] - self.startTime
        else:
            elapsedTime = self.endTime - self.startTime
        return elapsedTime.seconds


class OrderedDict(dict):
    """
    Custom dict subclass which stores each key that's added in an ordered
    list. However, it doesn't support initializing with a dict, which itself
    by definition will be unordered.  Provides basic ordered dict abilities
    on python older than 2.7
    """

    def __init__(self, *args, **kwargs):
        """Override __init__ to set up an orderedKeys list"""
        self.orderedKeys = []
        super(OrderedDict)


    def __setitem__(self, key, value):
        """Override __setitem__ to save key in the orderedKeys list"""
        if key not in self.orderedKeys:
            self.orderedKeys.append(key)
        super(OrderedDict, self).__setitem__(key, value)


    def items(self):
        """Override items to return ordered key,value pairs"""
        itemlist = []
        for key in self.orderedKeys:
            itemlist.append( (key, self.__getitem__(key)) )
        return itemlist


    def keys(self):
        """Override keys to return ordered key list"""
        return self.orderedKeys


    def values(self):
        """Override values to return values ordered by ordered key list"""
        valuelist = []
        for key in self.orderedKeys:
            valuelist.append( self.__getitem__(key) )
        return valuelist

