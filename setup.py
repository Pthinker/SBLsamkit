####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####

from distutils.core import setup
import SBLsamkit

NAME             = 'SBLsamkit'
VERSION          = SBLsamkit.__version__
DESCRIPTION      = 'SBLsamkit SAM Utilities'
AUTHOR           = 'Brian Sheehan'
AUTHOR_EMAIL     = 'sheehan@bgstech.com'
URL              = 'http://www.sbl-uk.org'
PLATFORM         = ['Unix', 'Linux', 'Mac OS X']
PACKAGES         = [ 'SBLsamkit', 
                   ]
SCRIPTS          = [ 'scripts/samkit.py',
                   ]
LONG_DESCRIPTION = """Python package which extends the capabilities of the
HTSeq package for filtering SAM streams based on multiple criteria"""

setup(name=NAME,
      version=VERSION,
      description=DESCRIPTION,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      url=URL,
      platforms=PLATFORM,
      long_description=LONG_DESCRIPTION,
      packages=PACKAGES,
      scripts=SCRIPTS,
     )
