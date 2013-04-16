####                       SBLsamkit module 4.19                            ####
#                                                                              #
#    SBLsamkit is part of a series of modules developed for Systems Biology    #
#    Laboratory UK, a community interest company involved in clinically        #
#    focussed research in immunology and molecular biology.                    #
#                                                                              #
####   Copyright 2010-2011 Systems Biology Laboratory, All Rights Reserved  ####



This directory contains some example samkit.ini files demonstrating various
ways you can define streams for a samkit.py run.  The .ini file *must* reside
in your working directory.  Each section in the samkit.ini file is started with
a section name.  The name can be any string, but avoid spaces or other 
non-alphanumeric characters.  Use the section name 'section01' if you wish
to pass --ini=samkit.ini without providing section name.  For example:

[set01]

You can then refer to that section when running samkit3.py with the --ini
argument, e.g.:

    samkit.py --if=your_input_file.sam --spn=200 --ini=samkit.ini,set01

You can define as many sections as you want in one or more samkit.ini files, 
provided each one has a unique name in any single file (different samkit.ini
files may share sections of the same name). Any stream defined within a section
must have a unique name, but you may use the same stream names in more than one
section (with the caveat that if you try to use both these sections in a single
run, an error will be generated (samkit does not allow multiple streams with
the same name):

[set01]
st1=all.fr,UUE,fr
st2=all.rf,UU,rf

[set02]
st1=uud_fr,UUD,fr
st2=uu_rf,UU,rf
