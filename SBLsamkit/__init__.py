####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

"""
samkit 4.x modules.  Designed based on previous SBL and SBLsam packages.

Copyright Systems Biology Laboratory (SBL) UK 2010-2013
All Rights Reserved

http://www.sbl-uk.org
"""

# Define module name and version:
__version_major__  = 4
__version_minor__  = 1
__version_update__ = 9
__version_bugfix__ = 0

if __version_bugfix__ == 0:
    __version__ = '%d.%d%s' % (__version_major__, 
                               __version_minor__,
                               str(__version_update__),
                              )
else:
    __version__ = '%d.%d%s.%s' % (__version_major__, 
                                  __version_minor__,
                                  str(__version_update__),
                                  str(__version_bugfix__),
                              )

from Processor import *

def version():
    print __version__

