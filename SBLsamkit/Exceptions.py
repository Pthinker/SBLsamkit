####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
Custom exceptions, each one accepts an optional second argument in the
initializer, a stream name.  If set, it will be incorporated into the error
message.

"""


class ParserException(Exception):
    """
    Base SBLsamkit exception which includes the stream name in the error
    message (the the stream is provided)

    """
    def __init__(self, message, stream=None):
        if stream is not None:
            message = "Stream %s: %s" % (stream, message)
        Exception.__init__(self, message)


class DuplicateStreamException(ParserException):
    """Custom exception for duplicate input stream names"""
    pass


class OptionParserException(ParserException):
    """Custom exception for parser errors"""
    pass


class LLParserException(ParserException):
    """Custom exception specifically for LL parser errors"""
    pass


class SSParserException(ParserException):
    """Custom exception specifically for SS parser errors"""
    pass


class iLRangeParserException(ParserException):
    """Custom exception specifically for iLrange parser errors"""
    pass


class IntRangeParserException(ParserException):
    """Custom exception specifically for NMrange,RPrange parser errors"""
    pass


class RnameListParserException(ParserException):
    """Custom exception specifically for RNAME list parser errors"""
    pass


class RnameListPairParserException(ParserException):
    """Custom exception specifically for RNAME pair list parser errors"""
    pass


class OperatorParserException(ParserException):
    """Custom exception specifically for operator parser errors"""
    pass

class OrderModeParserException(ParserException):
    """Custom exception specifically for order mode parser errors"""
    pass

class OrderMissingFormatParserException(ParserException):
    """Custom exception for error dues to missing output format mode"""
    pass

class PercentileParserException(ParserException):
    """Custom exception specifically for percentile argument parser errors"""
    pass

