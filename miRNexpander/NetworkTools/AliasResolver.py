#!/usr/bin/env python

from __future__ import print_function

#import MySQLdb
import os, sys
from miRNexpander.mWBBase import mWBBaseClass
from miRNexpander.DatabaseTools.DatabaseHandler import DatabaseHandler



class AliasResolver( mWBBaseClass ):
    """mediate translation between unique ids = shared ids = references and syonnyms = accessions = aliases"""

    def __init__( self, db_handler ):
        """initialize instance variables"""

        self._db = db_handler


############################################################
#### PUBLIC                                             ####
############################################################

    def setSpeciesRestriction( self, species_list ):
        """consider only the specified species for database queries"""

        if type( self._db.setAllowedSpecies( species_list ) ) == int:
            self._alert( "Failed to set species restriction, exiting." )
            sys.exit( 1 )
        else:
            self._spill( "Set species restriction to {}." . format( ", " . join( map( str, self._db._allowed_species ) ) ) )


    def getSpeciesRestriction( self ):
        """return the current restriction"""

        return self._db.getAllowedSpecies( )


    # - candidates: tuple, list, dict (keys will be used) or set
    def unalias( self, candidates, restrict = { "alias_types" : None }, silent = True, return_all = False ):
        """query database with aliases, return unique identifiers"""

        res = self._db.query_for_references( candidates, restrict = restrict, silent = silent )
        try:
            int( res )
            self._alert( "Error while calling embedded method." )
            return res
        except TypeError:
            pass

        if not silent:
            if res[ 1 ]:
                self._spill( "The following aliases remain orphaned:\n\t{}" . format( ", " . join ( sorted( res[ 1 ] ) ) ) )
            if res[ 2 ]:
                self._spill( "The following aliases were ambiguous:\n\t{}" . format( ", " . join( [  list( s )[ 0 ][ 1 ] for s in res[ 2 ]  ] ) ) )

        if return_all:
            return res
        else:
            return res[ 0 ]


    def aliases( self, candidates, restrictions = None, silent = True ):
        """query database with aliases for all aliases"""

        res = self._db.query_for_aliases( candidates, restrictions, silent = silent )
        try:
            int( res )
            self._alert( "Error while calling embedded method." )
            return res
        except TypeError:
            pass

        if not silent:
            if res[ 1 ]:
                self._spill( "The following identifiers did not produce any aliases under the given restrictions:\n\t", res[ 1 ] )

        return res[ 0 ]


    def annotate( self, candidates, silent = True ):
        """retrieve basic information about the queried molecules"""

        id_list = [ ]
        try:
            id_list += self.unalias( candidates[ "aliases" ], silent ).keys( )
        except KeyError:  # no key "aliases" - that is fine
            pass
        try:
            id_list += candidates[ "ids" ]
        except KeyError:  # no key "ids" - that is also fine
            pass

        if not id_list:
            return { }

        return self._db.query_for_annotations( id_list, silent )


############################################################
#### NOT SO PUBLIC                                      ####
############################################################

