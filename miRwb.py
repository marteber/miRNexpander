#!/usr/bin/python

from __future__ import print_function

from miRNexpander.DatabaseTools.DatabaseHandler import DatabaseHandler as DH
from miRNexpander.NetworkTools.NetworkCreator import NetworkCreator as NC
from miRNexpander.NetworkTools.CellDesignerIO import CellDesignerIO as CD
from miRNexpander.NetworkTools.InteractionLister import InteractionLister as IL
import argparse
import ConfigParser
import sys, os
import glob

###### COMMAND-LINE PARAMETER DEFINITION ######

# define command-line parameters ( aka arguments/options )
clp = argparse.ArgumentParser( description = "high-level cli to miRNexpander; provides all relevant functions" )

#clpg1 = clp.add_mutually_exclusive_group( )  # add group members with clpg.add_argument
#clpg1.add_argument( '-q', '--quiet', action = "store_true", help = "suppress all non-error output" )
#clpg1.add_argument( '-Q', '--mute', action = "store_true", help = "suppress all output" )

clp.add_argument( '-v', '--version', action = "version", version = "%(prog)s 1.1 (July 2014)" )
clp.add_argument( '-c', '--config', metavar = "CONFIG_FILE", help = "read configuration from this file", default = "setup.cfg" )

clp.add_argument( '-u', '--update', metavar = "DB", help = "try to update the local database; update all if no DB is specified", nargs = '*' )
clp.add_argument( '-w', '--write', metavar = "TARGET", help = "write updated sql statements to these files; write to standard locations if TARGET is empty", nargs = '*' )
clp.add_argument( '-r', '--read', metavar = "SOURCE_LIST", help = "read data from these files; one SOURCE_LIST per database, files must be separated with commas (,)", nargs = '+' )

clpg2 = clp.add_mutually_exclusive_group( )
clpg2.add_argument( '-a', '--alias', metavar = "MOL_ID | FILE", help = "print the reference ID(s) of the input MOL_ID(s)", nargs = '+' )
clpg2.add_argument( '-A', '--annotate', metavar = "MOLD_ID | FILE", help = "print the reference annotation of the input MOL_ID(s)", nargs = '+' )

clp.add_argument( '-l', '--flist', metavar = "FILE", help = "create an interaction list in Excel spreadsheet format", nargs = '+' )

clp.add_argument( '-p', '--celldesigner', metavar = "FILE", help = "parse the specified FILE into a network" )

clp.add_argument( '-b', '--build', metavar = "MOL_ID", help = "use these untyped seeds to build a network", nargs = '+' )
#clp.add_argument( '-B', '--files', metavar = "FILE", help = "read seeds from FILE(s); seeds must be separated with tabstops or newlines", nargs = '+' )
#clpg2.add_argument( '-o', '--online', metavar = "MIRNA", help = "like -b, but favor online requests for building", nargs = '+' )
#clpg2.add_argument( '-O', '--online-strict', metavar = "MIRNA", help = "like -o, but abort if an online request fails", nargs = '+' )
clp.add_argument( '-m', '--mirnas', metavar = "MIRNA_ID", help = "use these miRNA seeds to build a network", nargs = '+' )
#clp.add_argument( '-M', '--mirna-files', metavar = "FILE", help = "read miRNA seeds from these files", nargs = '+' )
clp.add_argument( '-n', '--name', metavar = "NETWORK_NAME", help = "use this name for the network (in filenames, e.g.)" )
clp.add_argument( '-d', '--depth', metavar = "DEPTH", help = "expand the network up to surrounding shell DEPTH", type = int, default = 1 )
clp.add_argument( '-s', '--species', metavar = "SPECIES", help = "restrict the network components to the specified species", nargs = '+' )

clp.add_argument( '-f', '--files', metavar = "FILE", help = "load seeds from these files and build a network", nargs = '+' )
clp.add_argument( '-g', '--filter', metavar = "FILE", help = "filter with gene Symbols from FILE", nargs = '+' )
#clpg3 = clp.add_mutually_exclusive_group( )
#clpg3.add_argument( '-i', '--include', metavar = "DB", help = "build only with the specified databases", default = "all" )
#clpg3.add_argument( '-e', '--exclude', metavar = "DB", help = "build with standard DBs excluding the specified databases", default = "none" )

nc = None  # pre-initialization of NetworkConstructor instance
il = None

###### SCRIPT PROPER ######

# parse parameters and react accordingly
parameters = clp.parse_args( )

# initialize configuration parser
conf = ConfigParser.ConfigParser( )
try:
    conf.readfp( open( parameters.config ) )
except IOError:
    sys.stderr.write( "Fatal: Could not read configuration file {!r}, quitting.\n" . format( parameters.config ) )
    sys.exit( 1 )

dh = DH( )
conn_error = dh.connect( conf )
if conn_error:
    sys.stderr.write( "Fatal: There's something wrong with the database connection: Error {:d}: {}\n" . format( conn_error.args[ 0 ], conn_error.args[ 1 ] ) )
    sys.exit( 1 )


### UPDATE (-u)

if parameters.update == [ ]:  # update goes first, empty list means update all
    #parameters.update = raw_input( "Really update all databases according to the configuration file? (y/N) " ) in ( "y", "yes" )
    #if parameters.update == True:
    #    print( "Trying to update all databases - please wait..." )
    dh.update_all( )
    parameters.update = False
elif parameters.update:
    print( "Trying to update database(s) {} - please wait..." . format( ", " . join( parameters.update ) ) )

if parameters.update:
    if parameters.write == [ ]:
        dh.update_db( parameters.update, True, parameters.read )
    elif parameters.write:
        dh.update_db( parameters.update, parameters.write, parameters.read, type( parameters.update ) == list )
    else:
        dh.update_db( parameters.update, source_list = parameters.read )


### ALIAS (-a) / Annotate (-A)

al_in = [ ]
if parameters.alias:
    al_in += parameters.alias
elif parameters.annotate:
    al_in += parameters.annotate

alias = set( )
if al_in:
    for f in al_in:
        try:
            stream = open( f )
            #parameters.unalias.extend( stream.read( ).splitlines( ) )
            alias |= {  l.strip( ).split( "\t" )[ 0 ] for l in stream  }
            #print( [ l.strip( ).split( "\t" )[ 0 ] for l in stream.read( ).splitlines( ) ] )
            stream.close( )
        except IOError:
            #print( "Not a file, using as string list.", file = sys.stderr )
            alias |= {  e for e in f.split( "," )  }

if alias:
    if not nc:
        nc = NC( dh )
    nc.set_species_restriction( parameters.species )
    ids = nc.unalias( alias, False )
    res = nc.annotate( {  "ids" : ids.keys( )  }, silent = False )

    if type( res ) == int:
        print( "An error occurred.", file = sys.stderr )
        sys.exit( res )

    if not res:
        print( "No valid reference IDs found.",  file = sys.stderr )
    if parameters.alias:
        for ref, alias in res.iteritems( ):
            print( alias[ "symbol" ] )
    else:
        for ref, alias in res.iteritems( ):
            print( "\t" . join( [  ids[ ref ]  ] + [  alias[ k ] for k in "symbol", "description"  ] ) )


###  FLIST (-l)
if parameters.flist:
    if not il:
        il = IL( dh )
    #print( il._db.getConfItem( "interdefs" ) )
    #print( il._db.getConfItem( "db_types" )[ "interactions" ] )
    #print( il.readIntseeds( parameters.flist ) )
    il.setSpeciesRestriction( parameters.species )
    if parameters.gfilter:
        gf = set( )
        for f in parameters.gfilter:
            try:
                gf |= set( open( f ).read( ).splitlines( ) )
            except IOError:
                print( "Unable to read file {!r}, skipping." . format( f ), file = sys.stderr )
    else:
        gf = None

    il.createExcel( *il.readIntseeds( parameters.flist ), tdir = "/media/sf_VM_Backup/Pia_validation", genefilter = gf, merge_cells = False )


### BUILD (-b/-B)

if ( parameters.build != None
    or parameters.mirnas != None
    #or parameters.proteins != None
    or parameters.files != None ):
    if parameters.files:
        for f in parameters.files:
            try:
                stream = open( f, "r" )
                for line in stream:
                    if line.strip( ) == "":
                        continue
                    try:
                        parameters.build += line.strip( ).split( "\t" )
                    except TypeError:
                        parameters.build = line.strip( ).split( "\t" )
                stream.close( )
            except IOError:
                sys.stderr.write( "Warning: File {!r} is not readable.\n" . format( f ) )
    parameters.build = list( set( parameters.build ) )  # remove duplicated entries
    seeds = { }
    seed_list = [ ]
    if parameters.mirnas != None:
        #seeds[ "miRNA" . lower( ) ] = parameters.mirnas
        seeds[ "unknown" ] = parameters.mirnas
        seed_list.extend( parameters.mirnas )
    if parameters.proteins != None:
        #seeds[ "protein" . lower( ) ] = parameters.proteins
        seeds[ "unknown" ] = parameters.proteins
        seed_list.extend( parameters.proteins )
    if parameters.build != None:
        seeds[ "unknown" ] = parameters.build
        seed_list.extend( parameters.build )
    #seeds = [ ]
    #if parameters.mirnas != None:
    #    seeds += parameters.mirnas
    #if parameters.proteins != None:
    #    seeds += parameters.proteins
    #if parameters.build != None:
    #    seeds += parameters.build

#if parameters.build:
    #print( "Building network\n\tup to shell {}\n\tfor {}\n\nPlease wait..." . format( parameters.depth, ", " . join( seeds ) ) )
    #print( "Building network\n\tup to shell {}\nPlease wait..." . format( parameters.depth ) )
    if not nc:
        nc = NC( dh )
    if not parameters.name:
        if len( seed_list ) == 1:
            parameters.name = seed_list[ 0 ]
        #if len( seeds ) == 1:
        #    parameters.name = seeds[ 0 ]
        else:
            parameters.name = "miRNA_centered_network"
    nc.setName( parameters.name )
    nc.set_species_restriction( parameters.species )
#    seeds = {  "unknown" : parameters.build  }
    nc._spill( "Building network\n\tup to shell {}\n\tfor {}\n\nPlease wait..." . format( parameters.depth, ", " . join( seed_list ) ) )
    if nc.createNetwork( seeds, parameters.depth ) == True:
        nc.writeSBML( )
        nc.writeXGMML( )
        nc.writeGML( )
        #nc.snapshot( "png" )
        #nc.writeCytoscape( )
    else:
        print( "No database entries for specified seeds list." )
#else:
#    print( "No network seeds specified, not building." )

### CELLDESIGNER (-p)
if parameters.celldesigner:
    if parameters.filter:
        gf = set( )
        for f in parameters.filter:
            try:
                gf |= set( open( f ).read( ).splitlines( ) )
            except IOError:
                print( "Unable to read file {!r}, skipping." . format( f ), file = sys.stderr )
    else:
        gf = None

    f = parameters.celldesigner
    cd = CD( dh )
    cd._spill( "Will try to convert CellDesigner to GML, please be patient.\n\t-> Parsing {}..." . format( f ) )
    cd.parse_file( f )
    manipulate = {
            "fold_to_symbol" : None,
            "link_subunits" : None,
            "bridge_reactions" : [ "TRANSCRIPTION", "TRANSLATION", "TRANSPORT" ],
            "bridge_nodes" : {
                "GENE" : None,
                "RNA" : None,
                "PHENOTYPE" : None,
                "ION" : None,
                "SIMPLE_MOLECULE" : None, #dict( graphics = dict( fill = [ "#CBFF66" . lower( ) ] ) ),
                "DEGRADED" : None,
            },
    }
    cd._spill( "\t-> Converting XML to graph..." )
    cd.xml_to_graph( manipulate )
    cd._spill( "\t-> Initial annotation from database..." )
    cd.annotate_graph( )
    cd._spill( "\t-> Writing output..." )
    cd.write_graph( f, ( "gml", ) )

sys.exit( 0 )


### draft area

