####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
Statistics functions designed to return lists of strings which can be easily
added to an InfoFile's ELN object
"""

import os
import math
import numpy
from Utils import *
from Operators import *


def SetPercentileValues(stream, opts, total):
    """
    This function analyzes a processed stream, specifically the array of
    insert lengths (sizes) in stream.stats['sizes'], and determines the
    percentile values in the passed-in dict.
    """
    try:
        Lmin,Lmax = stream.ops['iLrange']
    except IndexError:
        Lmin,Lmax = 0,numpy.inf

    for (key,value) in opts.pct.items():
        cutoff  = percentile(stream.stats['sizes'], value)
        opts.pct[key] = cutoff


def report_general_stats(opts, total):
    """
    Generate basic stream stats and return as an array of strings
    """

    # 4.14.4 - Set up list to track values needed for the "streams.csv" file
    streamSTATS = []
    streamSTATS.append(("Input", "n/a", total))

    output = ['','']
    output.append('Suffix %31s iLrange  Stream   Records  %% input:' % '')
    output.append('')
   
    output.append('%-51s %12d   %6.2f%%' % ('Input', total, 100.0))

    invalidCount = 0
    for stream in opts.orderedStreams:
        invalidCount += stream.stats['invalid']
        strSFX = '%s' % stream.streamSuffix
        # A blank suffix is allowed:
        if len(strSFX) > 0 and strSFX[0] == '.':
            strSFX = strSFX[1:]
        streamSTATS.append((strSFX, stream.name, stream.count))
        iLrange = iLrangeAsString(stream.ops['iLrange'])
        if total == 0:
            pct = 0.0
        else:
            pct = stream.count * 100.0 / total
        entry = '%-35s %10s  %-6s%10d   %6.2f%%' % (stream.streamSuffix,
                                                    iLrange,
                                                    stream.name,
                                                    stream.count,
                                                    pct)
        output.append(entry)

    # Invalid stream count (those which raised an exception on write)
    if total == 0:
        pct = 0.0
    else:
        pct = invalidCount * 100.0 / total
    entry = '%-35s %10s  %-6s%10d   %6.2f%%' % ('invalid',
                                                '0-inf',
                                                'n/a',
                                                invalidCount,
                                                pct)
    output.append(entry)

    return output, streamSTATS


def report_stream_summary(opts):
    outputs = ['',]
    outputs.extend(report_stream(None, opts))
    for stream in opts.orderedStreams:
        if stream.count > 0:
            outputs.extend(report_stream(stream, opts))
    return outputs


def report_stream(stream, opts):
    """
    Returns a list of strings describing LL and SS type distributions either
    from the global stats list of any stream if the first argument is None (the
    global stats will all be identical for any stream so we simply use the
    first stream in the opts.orderedStreams list) or from the specified stream 
    object.  In the latter case an additional column (% of input) is also
    computed and displayed.  The function layout is a little complicated, to
    allow two tables to be displayed side by side.
    """
    out = []

    hdr1 = 'LL           count:   %% of %-12s  UUE          count:   %% of %s'
    hdr2 = 'type:                  AA: %-12s  ss:                   UUE: %s'

    if stream is None:
        out.append('Input file summary:')
        out.append(opts.inputfile)
        out.append('')
        out.append(hdr1 % ('', ''))
        out.append(hdr2 % ('', ''))
    else:
        out.append('Output file summary:')
        # handle internal streams:
        if stream.name == 'n/a':
            out.append('%s      (%s)' % (stream.streamSuffix, 
                                         stream.name))
        else:
            out.append('%s%s    (%s)' % (stream.streamSuffix, 
                                         os.path.splitext(opts.inputfile)[1],
                                         stream.name))
        out.append('')
        out.append(hdr1 % ('  % of', '  % of'))
        out.append(hdr2 % ('input: ', 'input: '))
    out.append('')

    glstream = opts.orderedStreams[0]
    all = glstream.globalstats['AA']
    if stream is None:
        for tag in ('AA', 'UU', 'UUE', 'UUD', 'UX', 'UnM', 'UrM', 
                    'nMnM', 'rMrM', 'nMrM', 'nMA', 'rMA', 'XX', 
                    'fr', 'rf', 'ff', 'rr'):
            exec "%s = glstream.globalstats['%s']" % (tag, tag)
    else:
        for tag in ('AA', 'UU', 'UUE', 'UUD', 'UX', 'UnM', 'UrM', 
                    'nMnM', 'rMrM', 'nMrM', 'nMA', 'rMA', 'XX', 
                    'fr', 'rf', 'ff', 'rr'):
            exec "%s  = stream.stats['%s']" % (tag, tag)

    # Trap for divide by zero errors.
    try:
        pctUUE = 100.0/UUE				# Percent of UUE
    except ZeroDivisionError:
        pctUUE = 0.0			
    try:
        pctAA = 100.0/AA				# Percent of AA
    except ZeroDivisionError:
        pctAA = 0.0
    try:
        pctAll = 100.0/all				# Percent all records
    except ZeroDivisionError:
        pctAll = 0.0

    if stream is None:
        statsArray = (
          ['AA',  AA,  AA * pctAA,  '', 'aa', UUE, UUE * pctUUE],
          ['UU',  UU,  UU * pctAA,  '', 'fr', fr,  fr * pctUUE],
          ['UUE', UUE, UUE * pctAA, '', 'rf', rf,  rf * pctUUE],
          ['UUD', UUD, UUD * pctAA, '', 'ff', ff,  ff * pctUUE],
          ['UX',  UX,  UX * pctAA,  '', 'rr', rr,  rr * pctUUE],
        )
        for arr in statsArray:
            arr[1] = LocalizedInteger(arr[1])
            arr[5] = LocalizedInteger(arr[5])
            out.append('%-4s %14s %6.1f %13s %-3s %15s %6.1f' % tuple(arr))
    else:
        statsArray = (
          ['AA', AA, AA*pctAA, AA*pctAll, '', 'aa', UUE,UUE*pctUUE,UUE*pctAll],
          ['UU', UU, UU*pctAA, UU*pctAll, '', 'fr', fr, fr*pctUUE, fr*pctAll],
          ['UUE',UUE,UUE*pctAA,UUE*pctAll,'', 'rf', rf, rf*pctUUE, rf*pctAll],
          ['UUD',UUD,UUD*pctAA,UUD*pctAll,'', 'ff', ff, ff*pctUUE, ff*pctAll],
          ['UX', UX, UX*pctAA, UX*pctAll, '', 'rr', rr, rr*pctUUE, rr*pctAll],
        )
        for arr in statsArray:
            arr[1] = LocalizedInteger(arr[1])
            arr[6] = LocalizedInteger(arr[6])
            atp = tuple(arr)
            out.append('%-4s %14s %6.1f %6.1f %6s %-3s %15s %6.1f %6.1f' % atp)

    if stream is None:
        for tag in ('UnM', 'UrM', 'nMnM', 'rMrM', 'nMrM', 'nMA', 'rMA', 'XX'):
            count = glstream.globalstats[tag]
            counti = LocalizedInteger(count)
            pct = count * pctAA
            out.append('%-4s %14s %6.1f' % (tag, counti, pct))
    else:
        for tag in ('UnM', 'UrM', 'nMnM', 'rMrM', 'nMrM', 'nMA', 'rMA', 'XX'):
            count = stream.stats[tag]
            counti = LocalizedInteger(count)
            pct = count * pctAA
            globalpct = count * pctAll
            out.append('%-4s %14s %6.1f %6.1f' % (tag, counti, pct, globalpct))
    out.append('')
    return out


def report_rname_summary(opts):
    """Returns a list of strings describing the occurances of each RNAME
    within each stream
    """
    out = ['RNAME Statistics:', '']

    for stream in opts.orderedStreams:
        out.extend(report_rname_for_stream(stream, opts, 'RNsingle'))
        if stream.ops['LL'] == LL_UUD:
            out.extend(report_rname_for_stream(stream, opts, 'RNpairs'))
    return out


def report_rname_for_stream(stream, opts, key):
    """Returns a list os strings describing the occurance of all RNAMES for
    that stream.  There are two sets of stats, once for single RNAMES and
    one for RNAME pairs.  The second argument 'key' is the selector on which
    set of stats to use (dictionary key, may be 'RNsingle' or 'RNpairs')
    """

    total = 0
    for val in stream.stats[key].values():
        total += val
    if total == 0:
        return []

    out = []
    if key == 'RNsingle':
        # handle internal streams:
        if stream.name == 'n/a':
            out.append('%s      (%s) (%s)' % (stream.streamSuffix,
                                              stream.name,
                                              stream.ops['LL']))
        else:
            out.append('%s%s    (%s) (%s)' % (stream.streamSuffix,
                                              os.path.splitext(opts.inputfile)[1],
                                              stream.name,
                                              stream.ops['LL']))
        out.append('')
        if stream.ops['LL'] == LL_UUE:
            out.append('RNAME: %48s pair count: %% total:' % '')
        else:
            out.append('RNAME: %53s count: %% total:' % '')
    else:
        out.append('UUD pair:%46s pair count: %%top 10:' % '')
    out.append('')
    out.append('Total count: %41s %12d  100.0' % ('', total))
    keys = stream.stats[key].keys()
    keys.sort()

    if key == 'RNsingle':
        for rnameString in keys:
            pct = stream.stats[key][rnameString] * 100.0 / total
            out.append('%-54s %12d %6.1f' % (Ellipsis(rnameString),
                                             stream.stats[key][rnameString],
                                             pct))
    else:
        items = stream.stats[key].items()
        items.sort(cmp=lambda x,y: cmp(y[1], x[1]))
        topten = 0
        for i in range(0,min(len(items),10)):
            rnameString, value = items[i]
            topten += value

        for i in range(0,min(len(items),10)):
            rnameString, value = items[i]
            pct = value * 100.0 / topten
            rnamePair = rnameString.split(',')	# Split RNAME pair over 2 lines
            rp1 = Ellipsis(rnamePair[0]) + ','
            rp2 = Ellipsis(rnamePair[1])
            out.append(' %-53s %12d %6.1f' % (rp1, value, pct))
            out.append(' %s' % rp2)

    out.append('')
    return out


def report_size_distributions(opts):
    """
    Returns a list of strings describing size distribution stastics
    """
    # Array of output strings
    outputs = ['Size Distribution Tables',]

    # Report on each stream
    for stream in opts.orderedStreams:

        # Initialize
        total = 0
        cusum = 0.0
        desum = 100.0

        # Distribution limits
        iLmin,iLmax = stream.ops['iLrange']

        # Skip any stream with no size matches
        if len(stream.stats['sizes']) == 0:
            continue

        # Phase 1: sizes <= 1000bp
        if iLmin < 1000:
            oa, total, cusum, desum = genHisto(stream, 0, 1000, 40, '25bp', 
                                               '0-1K', total, cusum, desum)
            outputs.extend(oa)

        # Phase 2: sizes <= 20Kbp
        if iLmax > 1000 and total < stream.count:
            oa, total, cusum, desum = genHisto(stream, 1000, 20000, 19, '1Kbp', 
                                               '1K-20K', total, cusum, desum)
            outputs.extend(oa)

        # Phase 3: sizes <= 50Kbp
        if iLmax > 20000 and total < stream.count:
            oa, total, cusum, desum = genHisto(stream, 20000, 50000, 3, '10Kbp',
                                               '20K-50K', total, cusum, desum)
            outputs.extend(oa)

        # Phase 4: sizes > 50Kbp
        if iLmax > 50000 and total < stream.count:
            try:
                if max(stream.stats['sizes']) > 50000:
                    histmax = max(stream.stats['sizes']) + 2
                else:
                    histmax = 50002
            except ValueError:
                histmax = 50002
            oa, total, cusum, desum = genHisto(stream, 50000, histmax,1,'*all*',
                                               '50K+', total, cusum, desum)
            outputs.extend(oa)

    outputs.append('')
    return outputs


def genHisto(stream, minVal, maxVal, numBins, binsize, lengths, 
             total, cusum, desum):
    """Generates histogram table output for a given stream, minVal, maxVal,
    and number of bins.  binsize and lengths are output strings.  
    total, cusum and desum track total and percentage sums.
    """

    # numpy before 1.3 has a semantics warning (deprecation) for the
    # hisotgram() function and suggests adding the named argument "new=True" 
    # to the call. After 1.3 the deprecation was removed and later versions
    # indicate the "new-True" argument will be removed.  So we check on the
    # version and act appropriately
    npver = numpy.__version__.split('.')
    intver = 0
    try:
        if len(npver) > 0:
            intver += int(npver[0]) * 100
        if len(npver) > 1:
            intver += int(npver[1]) * 10
        if len(npver) > 2:
            intver += int(npver[2])
    except:
        intver = 0

    if 0 < intver < 130:
        histo,bins = numpy.histogram(stream.stats['sizes'],
                                     bins=numBins,
                                     range=(minVal,maxVal),
                                     new=True)
    else:
        histo,bins = numpy.histogram(stream.stats['sizes'],
                                     bins=numBins,
                                     range=(minVal,maxVal))
    oa = []
    oa.append('')
    oa.append('%s (%s): fragment lengths %s. Bin size %s' % (
              stream.streamSuffix, stream.name, lengths, binsize))
    oa.append('')
    oa.append('  Bin:           count:       %:    cusum:    desum:')
    oa.append('')

    # Should we examine iLmin and iLmax?
    try:
        check_MIN, check_MAX = stream.ops['iLrange']
    except IndexError:
        check_MIN, check_MAX = 0, maxVal + 1

    for i in range(0,len(histo)):
        num = int(round(histo[i]))
        bin = bins[i]
        try:
            pct = num * 100.0 / stream.count
        except ZeroDivisionError:
            pct = 0.0
        total += num
        cusum += pct
        desum -= pct
        cusum = min(cusum, 100.0)
        desum = max(desum, 0.0)
        if bin >= check_MIN and bin < check_MAX:
            if numBins == 1:
                oa.append('%5d+ %16d %7.2f%%  %7.2f%%  %7.2f%%' % (
                          minVal, num, pct, cusum, desum,))
            else:
                oa.append('%6d %16d %7.2f%%  %7.2f%%  %7.2f%%' % (
                          bin, num, pct, cusum, desum,))
    return oa, total, cusum, desum


def report_rname_stats(opts):
    """
    If one or more filters has an RNAME description, generate RNAME
    related stats for each relevant filter
    """
    output = ['']
    for stream in opts.orderedStreams:
        if len(stream.ops['rnames']) > 0:
            output.append('Matches for %s (%s) (%s):' % (stream.ops['suffix'],
                                                         stream.name,
                                                         stream.ops['LL']))
            output.append('')
            output.append('%-50s p_min,p_max:       count:' % 'RNAME')

            # RNAMES may contain an rname,p_min,p_max tuple for single RNAMEs,
            # or two such tuples, representing an range pair

            # RANGE pair
            if stream.ops['RNpairs']:
                for i in range(0, len(stream.ops['rnames']), 2):
                    rn0 = stream.ops['rnames'][i]
                    rn1 = stream.ops['rnames'][i+1]
                    p_min0 = rn0[1]
                    p_max0 = rn0[2]
                    if p_max0 == numpy.inf:
                        p_max0 = -1
                    p_min1 = rn1[1]
                    p_max1 = rn1[2]
                    if p_max1 == numpy.inf:
                        p_max1 = -1
                    rname = (rn0,rn1)
                    try:
                        count = stream.stats['rnames'][RnamePairAsKey(rname)]
                    except KeyError:
                        count = 0
                    output.append('%-42s %10s/%10s %12d'%(Ellipsis(rn0[0]),
                                                        '%d,%d' % (p_min0,p_max0),
                                                        '%d,%d' % (p_min1,p_max1),
                                                        count))

            # Single RNAME
            else:
                for rnitem in stream.ops['rnames']:
                    if type(rnitem) in (tuple,):
                        allrnames = rnitem
                    else:
                        allrnames = (rnitem,)
                    for rname in allrnames:
                        p_min = rname[1]
                        p_max = rname[2]
                        if p_max == numpy.inf:
                            p_max = -1
                        try:
                            count = stream.stats['rnames'][RnameAsKey(rname)]
                        except KeyError:
                            count = 0
                        output.append('%-52s %10s %12d'%(Ellipsis(rname[0]),
                                                         '%d,%d' % (p_min,p_max),
                                                           count))
            output.append('')
    return output


def compute_rname_stats(opts):
    """
    If one or more filters has an RNAME description, generate RNAME
    related stats for each relevant filter and return as a list of
    dicts containing individual rname/stream counts, along with an 
    ordered list of dict keys
    """

    # Construct lists rnames and suffixes
    rnameList     = UniqueRNAMEList(opts.orderedStreams)
    suffixList    = SuffixList(opts.orderedStreams, OP_CSV)

    # Construct a list RNAME match counts.  This takes the form of a dict
    # of dicts.  The top level key is taken from each of the RNAME,p_min,p_max 
    # tuples stored in rnameList, and toplevel values are themselves dicts,
    # where each key is taken from the suffixList, and each value is the
    # match count number.
    rnameStats = {}
    for rname in rnameList:
        # RNAME pair:
        if len(rname) == 2:
            rnameStats[RnamePairAsKey(rname)] = {}
            for suffix in suffixList:
                rnameStats[RnamePairAsKey(rname)][suffix] = 0
        # Single RNAME:
        else:
            rnameStats[RnameAsKey(rname)] = {}
            for suffix in suffixList:
                rnameStats[RnameAsKey(rname)][suffix] = 0

    # Set rname counts for each unique combination of RNAME, p_min, p_max and
    # suffix.  For RN pairs (rngp, or rngpf operators), the combination insead is
    # RNAME, p_min_1, p_max_1, RNAME, p_min_1, p_max_2
    for stream in opts.orderedStreams:
        suffix = stream.ops['suffix']
        if stream.ops['RNpairs']:
            for i in range(0,len(stream.ops['rnames']),2):
                rnitem0 = stream.ops['rnames'][i]
                rnitem1 = stream.ops['rnames'][i+1]
                rnitem = (rnitem0, rnitem1)
                try:
                    matchCount = stream.stats['rnames'][RnamePairAsKey(rnitem)]
                except KeyError:
                    matchCount = 0
                rnameStats[RnamePairAsKey(rnitem)][suffix] = matchCount
        else:
            for rnitem in stream.ops['rnames']:
                if type(rnitem) in (tuple,):
                    rnames = rnitem
                else:
                    rnames = (rnitem,)
                for rname in rnames:
                    try:
                        matchCount = stream.stats['rnames'][RnameAsKey(rname)]
                    except KeyError:
                        matchCount = 0
                    rnameStats[RnameAsKey(rname)][suffix] = matchCount

    return rnameStats


def percentile(data, percent):
    """
    Given an array of values and a percent value between 0 and 100,
    this function returns the percentile of the provided array.  For example, 
    for percent=50 this will return the median value

    This function is derived from a numerical recipe on the ActiveState web
    site, the link for whilch is available here:

    http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/

    We used this recipe so that we don't require all of scipy be installed,
    just for this one function.  Unlike the original recipe we sort the 
    provided array we get so it doesn't need to be sorted ahead of time

    Note that the value returned will always be rounded to the closest integer
    """

    sortedData = numpy.sort(data)

    key   = (len(sortedData)-1) * (percent/100.0)
    lower = math.floor(key)
    upper = math.ceil(key)

    try:
        if lower == upper: 
            return sortedData[int(key)]

        d0 = sortedData[int(lower)] * (upper - key)
        d1 = sortedData[int(upper)] * (key - lower)
        return int(round(d0 + d1))
    except IndexError:
        return 0

