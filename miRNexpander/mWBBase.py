#!/usr/bin/env python

from __future__ import print_function

import sys
from inspect import stack


class mWBBaseClass:
    """define some very basic useful functions"""

    def _spill( self, message ):
        """write message to STDERR"""
        print( message, file = sys.stderr )


    def _alert( self, message ):
        """add class reference to message, then write it to STDERR"""
        print( "{}: {}" . format( self.__class__.__name__, message ), file = sys.stderr )


    def _extalert( self, message ):
        """add class and function reference to message, then write it to STDERR"""
        print( "{}.{}( ): {}" . format( self.__class__.__name__, stack( )[1][3], message ), file = sys.stderr )

