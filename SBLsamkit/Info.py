####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
Standard Info file classes.  This sub-package handles creation, parsing and 
modification of the .info files
"""


import os
import re
import sys
import time
import types
import textwrap
from copy import deepcopy
from datetime import datetime

from Utils import LocalizedInteger
from SBLsamkit import __version__ as SBLversion


class Entry(object):
    """
    Define the basic ELN (Electronic Lab Notebook) entry list, which is
    a list of indented, text-wrapped strings bordered by a timestamp string

    """

    INDENT = '    '     # Standard indentation for ELN entries (4 spaces)
    MAXLINE = 80        # Maximum line length for ELN entries

    TS = re.compile("[0-9]{4}\.[0-9]{2}\.[0-9]{2} " +
                    "[0-9]{2}:[0-9]{2}:[0-9]{2}    ")

    def __init__(self, tstamp=None):
        """Initialize object with current timestamp is no timestamp provided"""
        if tstamp is None:
            tstamp = self.now
        self.tstamp = tstamp
        self.start()

    def start(self):
        """If the object doesn't already have a 'start' section, begin the 
        Entry with the timestamp and a blank line
        """
        if not hasattr(self, 'lines'):
            self.lines = ['%s    Start' % self.tstamp,]
            self.blank()

    def append(self, lines, unindented=True):
        """Appends the provided strings to the ELN entry.  The 'lines' argument
        can be a string or list/tuple of strings.  The optional unindented 
        argument defaults to True but can be set to false, which will 
        strip 4 a leading indent from each line before appending it.
        """
        if type(lines) in (types.TupleType, types.ListType):
            lines = list(lines)
        else:
            lines = [lines, ]
        for line in lines:
            self.lines.extend(self.__wrapline__(line, unindented))

    def prepend(self, lines, unindented=True):
        """Similar to the append() method, prepend inserts the proved string
        of string tuple/list to tbe beginning of the ELN entry, immediately
        after the timestamp and blank line set by the start() method
        """
        if type(lines) in (types.TupleType, types.ListType):
            lines = list(lines)
        else:
            lines = [lines, ]
        index = 2
        for line in lines:
            wrappedLines = self.__wrapline__(line, unindented)
            for wrapped in wrappedLines:
                self.lines.insert(index, wrapped)
                index += 1

    def blank(self):
        """Insert a blank line into the Entry"""
        self.lines.append('')

    def end(self):
        """If the object doesn't already finish with an 'end' section, append
        a blank line and the timestamp
        """
        line = '%s    End' % self.tstamp
        if line not in self.lines:
            self.blank()
            self.lines.append(line)

    def __repr__(self):
        """Combine the entry (a list of strings) into a single string with a 
        newline between each line, respecting correct indentation.  This allows
        the caller to simply print the object
        """
        
        out = []
        for line in self.lines:
            if line.strip() == '' and len(out) > 0 and out[-1].strip() == '':
                continue
            if self.TS.match(line):
                out.append(line)
            else:
                out.append('%s%s' % (self.INDENT, line))
        return '\n'.join(out)

    def __unindent__(self, line):
        """
        A replacement for the string class strip() method, which
        will remove no more than len(self.INDENT) number of spaces from the
        LHS of the incoming string. Any trailing spaces on the RHS will be 
        removed.
        """
        return re.sub('^ {0,%d}' % len(self.INDENT), '', line.rstrip())

    def __wrapline__(self, line, unindented=False):
        """
        Given a single line of text, wrap it into multiple lines by allowing
        no line to be longer than self.MAXLINE - len(self.INDENT) long, 
        preserving extra indentations if necessary.
        """
        if line.strip() == '':
            return ['', ]
        lines = []
        for wrappedline in textwrap.wrap(line, self.MAXLINE - len(self.INDENT)):
            if unindented:
                lines.append(wrappedline.rstrip())
            else:
                lines.append(self.__unindent__(wrappedline))
        return lines

    @property
    def now(self):
        """
        Returns a formatted string containing the current timestamp.
        """
        return datetime.now().strftime('%Y.%m.%d %H:%M:%S')


class RunEntry(Entry):
    """Subclass of the basic Entry class which adds suport for the program
       run subsection of the .info file ELN entry
    """

    def __init__(self, app, version, arglist, tstamp=None, entType='Automatic'):
        Entry.__init__(self, tstamp)
        self.app = app
        self.appver = version
        self.arglist = arglist
        self.entryType = entType

    def addRunSection(self, pct, pctLim, inputfile, outputdir, 
                      outputfiles, rangefiles, statsDict):
        """Sets up the run subsection and inserts it into the beginning of the
        ELN entry.  This allows us to set up the run subsection at any time
        before writing the file and it will still be at the beginning of the
        ELN entry.
        """
        self.pct         = pct
        self.pctLim      = pctLim
        self.inputfile   = inputfile
        self.outputfiles = outputfiles
        self.outputdir   = outputdir
        self.rangefiles  = rangefiles
        self.prepend(self.CommandLineSection +
                     self.InputFileSection +
                     self.OutputFileSection + 
                     self.OutputDirSection +
                     self.RangeFileSection +
                     self.PercentileSection + 
                     self.addRunStats(statsDict) +
                     self.ArgListSection)

    @property
    def ArgListSection(self):
        """Adds the argument list to the ELN entry"""
        out = []
        out.append('Arguments:')
        out.append('')
        out.extend(self.arglist)
        out.append('')
        return out

    @property
    def CommandLineSection(self):
        """Returns a list of strings which describe the command line used when
        the calling program was run
        """
        progdir,progname = os.path.split(sys.argv[0])
        out = []
        out.append('Program ID: %s VN: %s' % (self.app, self.appver))
        out.append('')
        out.append('Command line:')
        out.append('')
        out.append('%s%s %s' % (self.INDENT, progname, ' '.join(sys.argv[1:])))
        out.append('')
        return out

    @property
    def InputFileSection(self):
        """Returns the input file section as a list of strings"""
        input = self.inputfile
        out = []
        out.append("Input data files:")
        out.append('')
        out.append("%s%s  date modified: %s" % (self.INDENT, input, 
                                                self.filemodtime(input)))
        out.append('')
        return out

    @property
    def OutputFileSection(self):
        """Returns the output file section as a list of strings"""
        out = []
        out.append("Output data files:")
        out.append('')
        for file in self.outputfiles:
            out.append("%s%s" % (self.INDENT, file))
        out.append('')
        return out

    @property
    def OutputDirSection(self):
        """Returns the output directory section as a list of strings"""
        out = []
        out.append("Output directory: %s" % self.outputdir)
        out.append('')
        return out

    @property
    def RangeFileSection(self):
        """Returns any range files as a list of strings"""
        out = []
        if len(self.rangefiles) > 0:
            out.append("Range files:")
            out.append('')
            for (file, identifier) in self.rangefiles:
                out.append("%s%s  %s" % (self.INDENT, file, identifier))
            out.append('')
        return out
    

    @property
    def PercentileSection(self):
        """If the --percentile option was used, this property returns a list
        of strings describing the computed values of each computed variable
        """
        if self.pct is None or self.pctLim is None:
            return []
        if len(self.pct.keys()) == 0:
            return []
        hdr = 'Computed percentile values (%d records, limit %d bp):' % (
                                                   self.pctLim['sample_size'],
                                                   self.pctLim['iLcalcmax'])
        out = [ hdr, '' ]
        keys = self.pct.keys()
        keys.sort()
        for key in keys:
            value = self.pct[key]
            out.append('%s%s = %d' % (self.INDENT, key, value))
        out.append('')
        return out
            

    def addRunStats(self, statsDict):
        """Returns a list of strings containing the run statistics"""
        out = []

        # Total runtime is reported in seconds
        runtime = statsDict['runtime'] / 60.0
        out.append('Run time:                  %.1f minutes' % runtime)
        out.append('Run type:                  %s' % statsDict['runtype'])

        out.append('Warnings:                  %d' % len(statsDict['warnings']))
        for warning in statsDict['warnings']:
            out.append('%s%s' % (self.INDENT, warning))
        out.append('Error messages:            %d' % len(statsDict['errors']))
        for error in statsDict['errors']:
            out.append('%s%s' % (self.INDENT, error))
        totalrecords = LocalizedInteger(statsDict['count'])
        out.append('Input records processed:   %s record pairs' % totalrecords)
        out.append('Processors/cores used:     %d' % statsDict['numcpus'])
        out.append('')
        return out

    def filemodtime(self, filename):
        """Returns the file modification time as a string"""
        mtime = time.localtime(os.path.getmtime(filename))
        return time.strftime("%d %b %Y %H:%M:%S %Z", mtime)


class InfoFile(object):
    """
    Info file class, handles all aspects of the SBL .info spec 0.5.
    """
    # Tuple defining info file version
    INFOVERSION = (1,0)

    # ELN separator string
    SEPARATOR_STRING     = '-' * 40

    # Regular expressions to locate timestamps
    TS_START = re.compile("[0-9]{4}\.[0-9]{2}\.[0-9]{2} " +
                          "[0-9]{2}:[0-9]{2}:[0-9]{2}    Start$")
    TS_END   = re.compile("[0-9]{4}\.[0-9]{2}\.[0-9]{2} " +
                          "[0-9]{2}:[0-9]{2}:[0-9]{2}    End$")

    # End PR line
    ENDPR = ['@PR End']

    def __init__(self, app=None, version=None, filename=None, automatic=True):
        """Initialize new InfoFile object.
        If a filename is provided, attempt to parse it and incorporate it
        into the current object, otherwise create a new (empty) object
        """
        self.stats      = {}
        self.PRlines    = []
        self.ELNdict    = {}
        self.elnkeys    = []
        self.extraFiles = []
        self.configureApp(app, version)
        self.read(filename)
        self.addPR()

    def read(self, filename):
        if filename is None:
            return
        # Info file will not have the .gz extension
        if os.path.splitext(filename)[1] == '.gz':
            filename = os.path.splitext(filename)[0] 
        try:
            oldinfofile = filename + '.info'
            oldinfo = open(oldinfofile, 'r').readlines()
        except:
            return
        tstamp = None
        # Duplicate_Filter puts a leading line of dashes, so skip if found
        if oldinfo[0].startswith('------------------------'):
            oldinfo.pop(0)
        if oldinfo[0].startswith('@HD VN:'):
            for line in oldinfo:
                line = line.rstrip()
                if line.startswith('@PR') and line not in self.ENDPR:
                    self.PRlines.append(line)
                elif self.TS_START.match(line):
                    tstamp = line.replace('    Start', '')
                    if tstamp not in self.elnkeys:
                        self.elnkeys.append(tstamp)
                        self.ELNdict[tstamp] = Entry(tstamp)
                elif self.TS_END.match(line):
                    self.ELNdict[tstamp].end()
                    tstamp = None
                elif tstamp is not None:
                    self.ELNdict[tstamp].append(line, unindented=False)
        else:
            tstamp = datetime.fromtimestamp(os.stat(oldinfofile).st_mtime)
            tstamp = tstamp.strftime('%Y.%m.%d %H:%M:%S')
            # Non-compliant info files assume first two strings on first line
            # define progname and version.
            try:
                progname = oldinfo[0].split()[0]
                progvers = oldinfo[0].split()[1].replace('(','').replace(')','')
            except IndexError:
                progname = 'unknown'
                progvers = 'unknown'
            self.PRlines.append(self.customPR(tstamp, progname, progvers))
            if tstamp not in self.elnkeys:
                self.elnkeys.append(tstamp)
                self.ELNdict[tstamp] = Entry(tstamp)
            for line in oldinfo:
                line = line.rstrip()
                self.ELNdict[tstamp].append(line, unindented=False)
            self.ELNdict[tstamp].end()
            tstamp = None

    def customPR(self, timestamp, name='', version=''):
        """Generates a standard @PR line from a basic timestamp and
        program name and version.  If either string is empty, set default values
        of "unknown".  Returns formatted string.
        """
        if name == '': name = 'unknown'
        if version == '': version = 'unknown'
        return '@PR %s  %s %s' % (timestamp, name, version)

    def configureApp(self, app, version):
        """Configures application name and version from sys.argv[0] if no
        app and/or version were provided to the constructor
        """
        if app is None:
            appfile = os.path.split(sys.argv[0])[1]
            app = os.path.splitext(appfile)[0]
        if version is None:
            version = SBLversion
        self.app = app
        self.appver = version

    def configureHeader(self, outputdir, outputfiles, recordCount):
        """Configures values for the @DF and @RN header lines"""
        # Extract sorted list of unique output filenames
        self.streamfiles = list(set(outputfiles))
        self.streamfiles.sort()
        self.count = recordCount

    def setInput(self, filenames):
        if type(filenames) in (types.TupleType, types.ListType):
            self.inputfiles = tuple(filenames)
        else:
            self.inputfiles = (filenames,)
        self.ELN.addInputFiles(self.inputfiles)

    def setOutput(self, directory, filenames):
        self.outputdir = directory
        if type(filenames) in (types.TupleType, types.ListType):
            self.outputfiles = tuple(filenames)
        else:
            self.outputfiles = (filenames,)
        self.ELN.addOutputFiles(self.outputfiles)
        self.ELN.addOutputDir(self.outputdir)

    def addPR(self, string=None):
        if string is None:
            string = '%s %s' % (self.app, self.appver)
        self.PRlines.append('@PR %s  %s' % (self.now, string))

    def append2ELN(self, lines):
        if type(lines) is types.StringType:
            lines = lines.split('\n')
        self.ELN.append(lines)

    def appendArgList(self):
        self.ELN.addArgList()

    def configureRunSection(self, opts, outputdir, outfiles, rangefiles):
        pct    = opts.pct
        if hasattr(opts, 'pctLim'):
            pctLim = opts.pctLim
        else:
            pctLim = None
        infile = opts.inputfile
        self.stats['numcpus'] = opts.nproc
        self.ELN.addRunSection(pct, pctLim, infile, outputdir, outfiles, 
                               rangefiles, 
                               self.stats)

    def configureExtraFilesSection(self, fileList):
        fcontent = []
        allExtra = tuple(self.extraFiles) + tuple(fileList)
        if len(allExtra) > 0:
            fcontent.append('')
            fcontent.append('')
            fcontent.append('Additional Info Files:')
            fcontent.append('')

        for file,identifier in allExtra:
            ifile = '%s.info' % file
            if os.path.exists(ifile):
                fcontent.append('')
                fcontent.append('File: %s; %s' % (file, identifier))
                fcontent.append('Begin >>')
                fcontent.extend(open(ifile).readlines())
                fcontent.append('<< End')
                fcontent.append('')

        self.append2ELN(fcontent)

    def openELN(self, arglist=[], automatic=True):
        self.ELNdict[self.now] = RunEntry(self.app, self.appver, arglist,
                                          self.now)

    def closeELN(self):
        self.ELN.end()

    def write(self, filename=None):
        self.closeELN()
        if filename is not None:
            fp = open(filename, 'w')
        else:
            fp = sys.stdout
        print >> fp, self.HD
        print >> fp, self.DF
        print >> fp, self.RN
        print >> fp, self.PR
        print >> fp, self.SEP
        print >> fp, self.ELN
        print >> fp, self.SEP
        for key in self.elnkeys:
            print >> fp, self.ELNdict[key]
            print >> fp, self.SEP
        if filename is not None:
           fp.close()


    @property
    def ELN(self):
        if self.now not in self.ELNdict.keys():
            self.openELN()
        return self.ELNdict[self.now]

    @property
    def HD(self):
        return "@HD VN: %s%8s%s %s" % (self.InfoVersion,
                                       ' ', self.app, self.appver)

    @property
    def DF(self):
        dfline = '@DF ' + ' '.join(self.streamfiles)
        return dfline.rstrip()

    @property
    def RN(self):
        return '@RN %d records    %d record pairs' % (self.count*2,self.count)

    @property
    def PR(self):
        return '\n'.join(self.PRlines + self.ENDPR)
        
    @property
    def SEP(self):
        return self.SEPARATOR_STRING

    @property
    def InfoVersion(self):
        return '%d.%d' % self.INFOVERSION

    @property
    def now(self):
        """
        Returns a formatted string containing the current timestamp.
        """
        if not hasattr(self, 'tstamp_now'):
            self.tstamp_now = datetime.now().strftime('%Y.%m.%d %H:%M:%S')
        return self.tstamp_now

    def copy(self):
        """
        Returns a deep copy of itself.
        """
        return deepcopy(self)

