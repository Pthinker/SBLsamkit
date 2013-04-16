####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""Command line and INI file parser classes and functions"""

import re
import sys
from optparse import OptionParser as BasicOptionParser
from ConfigParser import SafeConfigParser
from SBLsamkit import __version__
from Utils import Error, Tokenize, IterateInPairs, Integer, OrderedDict
from Operators import *
from Stream import SAMStream
from Exceptions import *
from multiprocessing import cpu_count


def parse_int(option, opt_str, value, parser):
    """
    Integer converter callback which uses Utils.Integer() to convert, e.g.
    100K to 100000
    """
    orig_value = value
    try:
        value = Integer(value)
        exec "parser.values.%s = %s" % (option.dest, value)
    except ValueError:
        raise OptionParserException("Invalid int value '%s' for '%s' option" % 
                                    (orig_value, option))

def parse_nproc(option, opt_str, value, parser):
    """
    Enhanced integer converter callback which accepts integer values and
    then compares them against a valid range (minimum 1, maximum being the
    total number of CPUs in the system minus 2).
    """
    ncpus = cpu_count()
    orig_value = value
    try:
        value = int(value)
        if value < 1 or value > ncpus:
            raise OptionParserException("Invalid np '%d': Must be an integer between 1 and %d (or the word 'all')" % (value, ncpus))

        exec "parser.values.%s = %d" % (option.dest, value)
    except ValueError:
        if value == 'all':
            value = ncpus
            exec "parser.values.%s = %d" % (option.dest, value)
        else:
            raise OptionParserException("Invalid np '%s': Must be an integer between 1 and %d (or the word 'all')" % (orig_value, ncpus))
    except ValueError:
            raise OptionParserException("Invalid np '%s': Must be an integer between 1 and %d (or the word 'all')" % (orig_value, ncpus))


def parse_stream(option, opt_str, value, parser):
    """
    Main stream parser function, instantiates a stream class based on the 
    requested stream.  In addition it stores the stream parser objects in 
    order in the parser.values.orderedStreams attribute to allow streams to 
    be processed in order.

    Streams may defined in the following way (operators are optional):
    --stN=suffix,[LL,][ss,][NMf,NMrange,][RPf,RPrange,][phredf,min,max,][iL,min,max,][rng,(range list) | rnpg(range pair list)]|[rngf,rangefile]|[rngfprangepairfile][,operator,operator...]
    """

    # Instantiate the appropriate stream object and append to the list of
    # ordered streams (processing order must match the order the streams were
    # defined, so you'll often see constructs like this:
    #
    #    for stream in opts.orderedStreams:
    #         ... do something()
    #
    stream = SAMStream(option.dest, value)
    if hasattr(parser.values, 'orderedStreams'):
        parser.values.orderedStreams.append(stream)
    else:
        parser.values.orderedStreams = [stream,]


class INIParser(SafeConfigParser):
    """
    An extended form of SafeConfigParser which stores the order of options for
    each section
    """

    ALLOWED_STREAMS = ('st',)

    def __init__(self, inifile, defaults=None):
        SafeConfigParser.__init__(self, defaults)
        self.inifile = inifile
        # Make the _sections list an ordered dict, so that the section
        # names occur in order.
        self._sections = OrderedDict()
        self.orderedKeys = {}

    def readfp(self, filehandle=None):
        if filehandle is None:
            filehandle = open(self.inifile,'r')
        SafeConfigParser.readfp(self, filehandle)
        for section in self._sections.keys():
            for key,value in self._sections[section].items():
                if key and value:
                    self._sections[section][key] = value.replace('\n', '')

    def optionxform(self, optionstr):
        """
        Overrides SafeConfigParser.optionxform() method to prevent converting
        option names to lowercase, and to store defined keys in an ordered list.
        This works because the SafeConfigParser.readfp() method internally
        calls SafeConfigParser._read() to parse the .ini file, and each new
        option is passed to optionxform() to process (by default converting it
        to lowercase).  So each time an option string is passed to optionxform,
        it happens in the order those options are defined in the .ini file.
        """
        # Store the current stream in the correct ordered section 
        currentSection = self._sections.keys()[-1]
        try:
            self.orderedKeys[currentSection].append(optionstr)
        except KeyError:
            self.orderedKeys[currentSection] = [optionstr,]
        return optionstr


    def items_as_arguments(self, section):
        """
        This method returns the ordered list of options as a list of properly
        formatted strings as the user would enter on the command-line
        """
        arguments = []
        for orderedKey in self.orderedKeys[section]:
            for key,value in self.items(section):
                if key == orderedKey:
                    generic_key = re.sub('[0-9]*$', '', key)
                    if generic_key in self.ALLOWED_STREAMS:
                        arguments.append('--%s=%s' % (key, value))
        return arguments


class OptionParser(BasicOptionParser):
    """
    This class adds ARGUMENT_LIMIT options for a given option 
    name.  So for example if you add an option --foo, the option parser 
    will not support --foo directly, but WILL support --fooN where N is any
    number in the range 1 to ARGUMENT_LIMIT, inclusive
    """

    ARGUMENT_LIMIT = 99

    def __init__(self, **kwargs):
        """
        Call the parent class initializer, then also initialize a list
        of "generic forms" of each argument to be added, which will be
        used to trim down the help output.
        """
        BasicOptionParser.__init__(self, **kwargs)
        self.generic_arguments = []
        self.iniFiles = []


    def add_option(self, *args, **kwargs):
        """
        Store the argument in our generic_arguments list, before handing 
        off to the parent add_option method.  This allows us to record the
        option name as one which will be output by the --help option: see
        self.format_option_help
        """
        if 'dest' in kwargs.keys():
            dest = kwargs['dest']
            if dest not in self.generic_arguments:
                self.generic_arguments.append(args[0])
        BasicOptionParser.add_option(self, *args, **kwargs)


    def add_multi_option(self, *args, **kwargs):
        """
        Add the option as normal using the parent class add_option method.
        In addition, add a numbered set of options corresponding to the 
        original option with an integer appended to it.
        """
        dest = kwargs.pop('dest')
        if dest not in self.generic_arguments:
            self.generic_arguments.append("--%s" % args[0])
        for i in range(1,self.ARGUMENT_LIMIT + 1):
            kwargs['dest'] = '%s%d' % (dest,i)
            args = ('--%s%d' % (dest,i),)
            self.add_option(*args, **kwargs)


    def format_option_help(self, formatter=None):
        """
        Only report on the generic form of an argument, so to avoid
        repetition of the help string for a given argument set configured
        self.add_multi_option.
        """
        # Go through the options list, adding the generic form of each to a
        # a list which will be consulted when displaying help output
        self.generic_arguments.append('--help')
        self.generic_arguments.append('--version')
        trimmed_options = []
        for option in self.option_list:
            optstring = option.get_opt_string()
            # Store the automatic --help and --version options
            if optstring == '--version' or optstring == '--help':
                trimmed_options.append(option)
                continue
            # Construct generic form of the argument name
            generic_form = re.sub('[0-9]+$', '', optstring)
            if generic_form in self.generic_arguments:
                option.dest = re.sub('[0-9]+$', '', option.dest)
                # Store a generic form of the argument string, substituting an 
                # N for an any trailing integer value, e.g. --ss1 -> --ssN
                for i in range(0,len(option._long_opts)):
                    long_opt = option._long_opts[i]
                    option._long_opts[i] = re.sub('[0-9]+$', 'N', long_opt)
                trimmed_options.append(option)
                self.generic_arguments.remove(generic_form)
        self.option_list = trimmed_options

        # Run parent class format_option_help method
        help_string = BasicOptionParser.format_option_help(self, formatter)

        # Now split the array into lines and locate the first option describing
        # a stream.  Insert a label before that line, and append additional help
        # info the end of the array, then merge the array of lines back into a 
        # string and return it.
        help_array  = help_string.split('\n')
        str_label = "Streams, N in range 1-%d: " % self.ARGUMENT_LIMIT
        for i in range(0, len(help_array)):
            if re.search('--stN=', help_array[i]):
                help_array.insert(i, '\n' + str_label)
                break

        # Strip redundant arguments from help string.  For example, the
        # help string --stN=ST becomes --st=
        # Store the index to the "Streams" help info as we'll be inserting
        # additional help info *before* this point
        firstI = -1
        streamStrings = []
        for i in range(0, len(help_array)):
            if re.search('stN=ST', help_array[i]):
                help_array[i] = re.sub('stN=ST\s*', 
                                       'stN=              ', 
                                       help_array[i])
                if firstI < 0:
                    firstI = i-1

        # Insert operator help strings into correct location
        operatorlist = list(OPERATORHELP)
        operatorlist.reverse()
        for help in operatorlist:
            help_array.insert(firstI, help)
        help_array.insert(firstI, '')

        # Add final blank line and return
        help_array.append('')
        return '\n'.join(help_array)


    def parse_args(self, args=None, values=None, verbose=True):
        """
        This method overrides the base class parse_args by performing an
        initial iteration through the input arguments (defaults to sys.argv[1:]
        in the same fashion as in the parent class), and storing them in
        an ordered list.  This gets returned the the caller as the third
        element of a tuple, the other two being the usual opts and args returned
        by the parent class parse_args() method.  This method also "pre-scans"
        the options list for a "--ini" flag, and parses the '*.ini'
        file for a set of command-line arguments to insert into sys.argv.
        This allows command-line options to be stored in a .ini file.  Assuming
        a vaid section is located in a *.ini file, the extra arguments
        which will be prepended to the existing list of command line args,
        is also returned, as the fourth and final item in the returned tuple.
        """
        ordered_args = []
        sectionargs = []
        sections = []
        if args is None:
            args = sys.argv[1:]

        # Pre-scan argument list for an --ini argument
        for arg_as_input in args:
            if arg_as_input.startswith('--ini'):
                i_opt, i_arg = BasicOptionParser.parse_args(self, 
                                                            [arg_as_input,],
                                                            values)
                sections.extend(i_opt.ini)
        if len(sections) > 0:
            sectionargs = self.__parse_ini_files__(sections)

        # Now re-scan arguments, including any additional arguments which
        # were inserted from a valid section in the .ini file
        args += sectionargs
        for arg_as_input in args:
            if not arg_as_input.startswith('--'):
                continue
            arg_as_input = arg_as_input.replace('--', '')
            ordered_args.append(re.sub('=.*$', '', arg_as_input))
        p_opts, p_args = BasicOptionParser.parse_args(self, args, values)

        # Verify a stream was defined
        if not hasattr(p_opts, 'orderedStreams'):
            Error('You must define at least 1 filter stream (e.g. --st1=...)')

        # Display summary output
        if verbose:
            print
            print "Requested and implied arguments:"
            for arg in sys.argv[1:] + sectionargs:
                print '  %s' % arg
            print

        return p_opts, p_args, ordered_args, sectionargs, self.iniFiles


    def __parse_ini_files__(self, sections):
        new_args = []
        defaultSection = 'section01'
        for item in sections:
            itemelements = item.split(',')
            if len(itemelements) == 1:
                filename,section = itemelements[0],defaultSection
            elif len(itemelements) == 2:
                filename,section = itemelements[0],itemelements[1].lower()
            else:
                Error("--ini=%s is an invalid option! %s" % item)
            try:
                config = INIParser(filename)
                config.readfp()
                self.iniFiles.append(filename)
            except IOError:
                Error("--ini=%s requested but no %s file found" % (item,
                                                                   filename))
            if section not in config.sections():
                Error("No such section '%s' found in %s" % (section,
                                                            filename))

            for option in config.items_as_arguments(section):
                new_args.append(option)
        return new_args
        

