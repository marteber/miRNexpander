#!/usr/bin/env python

from __future__ import print_function

from miRNexpander.mWBBase import mWBBaseClass
import MySQLdb
#from webcolors import name_to_hex, rgb_to_hex
import ConfigParser
import sys, os
from inspect import stack


class DatabaseConnector( mWBBaseClass ):
    """handle connections and queries to the database"""


    def __init__( self ):
        """initialize values"""
        self._config = None
        self._cursor = None
        self._serv = None

        self._db_keys = set( [ "INDEX", "UNIQUE" ] )

        self._cytoscape_shapes = [  "DIAMOND", "ELLIPSE", "HEXAGON",
                                    "OCTAGON", "PARALLELOGRAM", "RECTANGLE",
                                    "RHOMBUS", "ROUNDED_RECTANGLE", "TRIANGLE", "VEE"
                                    ]


    def connect( self, conf ):
        "establish the database connection"""

        self._config = conf
        self._conf_dict = self._read_conf( self._config )

        host = self._config.get( 'database', 'host' )
        user = self._config.get( 'database', 'user' )
        passwd = self._config.get( 'database', 'passwd' )
        db = self._config.get( 'database', 'database' )

        try:
            self._serv = MySQLdb.connect( host, user, passwd, db )
        except MySQLdb.Error as e:
            return e
        self._cursor = self._serv.cursor( )
        #self._cursor.connection.autocommit( True )


    def _read_conf( self, conf ):
        """read the configuration for the list of specified tables into a dictionary"""

        def tc( t ):  # type correction hack (type hierarchy not implemented yet)
            #if t in ( "mirna", "pri-mirna" ):
            #    return t
            #else:
            #    return "protein"
            return t

        _dbprefix = "db_"
        result = {
                "suffixes" : {  "format" : "_format", "tmp" : "_tmp", "org" : "_db"  },
                "query_dbs" : set( ),
                "interdefs" : dict( ),
                "inter_fn" : dict( ),
                }
        _fmt = result[ "suffixes" ][ "format" ]
        _tmp = result[ "suffixes" ][ "tmp" ]

        _fn_template = "{rdate}_{filtered}{searched}_{searchedtag}_of_{ldate}_listed_{listed}_{species}"

        try:
            result.update( conf.items( "miRNexpander" ) )
            #result[ "known_dbs" ] = dict( [ ( t.lower( ), t ) for t in result[ "known_dbs" ].split( ) ] )
            result[ "known_dbs" ] = dict( [ ( t[ 3: ].lower( ), t[ 3: ] ) for t in conf.sections( ) if t[ :3 ] == _dbprefix ] )
            result[ "moltypes" ] = { k.lower( ) for k in result[ "moltypes" ].split( ) }
        except ( ConfigParser.NoSectionError, KeyError ):
            self._extalert( "Entries for section {} incomplete." . format( "mirnaworkbench or table_setup" ) )
            return 1

        # archetypes deciphering !TODO!
        helper = dict( conf.items( "archetypes" ) )
        result[ "moltypes" ] |= { k.lower( ) for k in helper.keys( ) + helper.values( ) }  # temporary

        # interaction definitions deciphering
        d = {  # necessary to keep format keys
                "rdate" : "{rdate}",
                "ldate" : "{ldate}",
                "species" : "{species}",
                "filtered" : "{filtered}",
            }
        for tup in conf.items( "interdefs" ):
            d.update( zip( ( "src", "tgt" ), tup[ 1 ].split( ) ) )
            result[ "interdefs" ][ tup[ 0 ] ] = {
                                                    "src" : { tc( k.lower( ) ) for k in d[ "src" ].split( ":" ) },
                                                    "tgt" : { tc( k.lower( ) ) for k in d[ "tgt" ].split( ":" ) },
                                                }
            for role, oppos, tag in ( "src", "tgt", "targets" ), ( "tgt", "src", "regulators" ):
                d.update( {
                            "listed" : "+" . join( map( lambda x: x + "s", d[ role ].split( ":" ) ) ),
                            "searched" : "+" . join( d[ oppos ].split( ":" ) ),
                            "searchedtag" : tag,
                            } )
                result[ "inter_fn" ][ tup[ 0 ] + "." + role ] = _fn_template . format( **d )
        #print( result[ "inter_fn" ] ); sys.exit( )

        result[ "tables" ] = dict( [ ( t, {  "UNIQUE": [ ], "INDEX": [ ]  } ) for t in result[ "known_dbs" ] ] )
        result[ "db_types" ] = {
                "annotations" : dict(  ),
                "aliases" : dict( ),
                "interactions" : dict( [ ( k, dict( ) ) for k in result[ "moltypes" ] ] ),
                }
        result[ "db_mapper" ] = dict( )
        result[ "sql_mapper" ] = dict( )
        sql = result[ "sql_mapper" ]  # alias into dictionary for easier access
        dbt = result[ "db_types" ]

        for tl, table in result[ "known_dbs" ].items( ):
            tab_id = _dbprefix + table
            t = result[ "tables" ][ tl ]  # alias into dictionary for easier access
            result[ "db_mapper" ][ tab_id ] = tl

            try:
                t.update( conf.items( tab_id ) )  # update dict with values from the corresponding section
                t[ "id" ] = tab_id
                try:
                    t[ "update" ] = int( t[ "update" ] )
                except ValueError:
                    if t[ "update" ] not in ( "yes", "y" ):
                        t[ "update" ] = 0
                t[ "paths" ] = [ os.path.join( t[ "dir" ], k ) for k in t[ "sources" ].split( ) ]
                #t[ "write" ] = t[ "write" ].split( )[ 0 ]  # TODO: rectify list handling in derived class
                t[ "write" ] = "txt"

                if "tablename" not in t or t[ "tablename" ] == "":
                    t[ "tablename" ] = table
                t[ "name" ] = t[ "tablename" ]
                sql[ tl ] = t[ "name" ]

                t[ "columns" ] = t[ "columns" ].split( )
                # interpret type of database based on parameters "src" and "tgt", "alias", or "annotation"
                if set( ( "src", "tgt" ) ) <= set( t ):  # we found an interaction database
                    t[ "src" ], t[ "tgt" ] = tc( t[ "src" ].lower( ) ), tc( t[ "tgt" ].lower( ) )
                    if set( ( t[ "src" ], t[ "tgt" ] ) ) <= result[ "moltypes" ]:
                        if t[ "src" ] == t[ "tgt" ]:
                            try:
                                dbt[ "interactions" ][ t[ "src" ] ][ "two-way" ].append( t[ "name" ] )
                            except KeyError:
                                dbt[ "interactions" ][ t[ "src" ] ][ "two-way" ] = [ t[ "name" ] ]
                        else:
                            try:
                                dbt[ "interactions" ][ t[ "src" ] ][ "src" ].append( t[ "name" ] )
                            except KeyError:
                                dbt[ "interactions" ][ t[ "src" ] ][ "src" ] = [ t[ "name" ] ]
                            try:
                                dbt[ "interactions" ][ t[ "tgt" ] ][ "tgt" ].append( t[ "name" ] )
                            except KeyError:
                                dbt[ "interactions" ][ t[ "tgt" ] ][ "tgt" ] = [ t[ "name" ] ]
                    else:
                        self._alert( "Warning: Config section {!r} has undefined types in {!r} and/or {!r}. Typo?" . format( t[ "id" ], "src", "tgt" ) )
                        sys.exit( 1 )
                    if t[ "retrieve" ] == "-":
                        t[ "retrieve" ] = [ ]
                    elif t[ "retrieve" ]:
                        t[ "retrieve" ] = [ int( k ) for k in t[ "retrieve" ].split( ) ]
                    else:  # default for interaction databases is columns 3 and 4
                        t[ "retrieve" ] = [ 3, 4 ]
                elif set( ( "annotation", "retrieve" ) ) <= set( t ):
                    dbt[ "annotations" ][ t[ "annotation" ].lower( ) ] = t[ "name" ]
                    if t[ "retrieve" ] == "-":
                        t[ "retrieve" ] = [ ]
                    elif t[ "retrieve" ]:
                        t[ "retrieve" ] = [ int( k ) for k in t[ "retrieve" ].split( ) ]
                    else:  # default for annotation databases is all columns
                        t[ "retrieve" ] = [ k + 1 for k in xrange( len( t[ "columns" ] ) ) ]
                    #new_name = t[ "name" ] + "_xrefs"
                    #add_table[ new_name ] = {
                    #        "alias" : None,
                    #        "columns" : [ new_name + "001:ref_key:I:n", new_name + "002:ref_type:40:n", new_name + "003:Alias:100:ni" ],
                    #        }
                elif set( ( "alias", ) ) <= set( t ):
                    dbt[ "aliases" ][ t[ "alias" ].lower( ) ] = t[ "name" ]

                # interpret column definitions and set up dependency resolution framework
                #try:
                #    t[ "depends" ] = t[ "depends" ].split( )
                #except KeyError:
                #    #t[ "depends" ] = ( " " . join( [ "db_gene_aliases" ] * 2 ) ).split( )
                #    t[ "depends" ] = ( " " . join( [ "db_gene_xrefs" ] * 2 ) ).split( )
                for i in xrange( len( t[ "columns" ] ) ):
                    col = t[ "columns" ][ i ]
                    ( c_id, colname, coltype, colprops ) = col.split( ':' )
                    try:
                        t[ "col_ids" ].append( c_id )
                        t[ "col_ids_with_tmps" ].append( c_id )
                    except KeyError:
                        t[ "col_ids" ] = [ c_id ]
                        t[ "col_ids_with_tmps" ] = [ c_id ]
                    try:
                        t[ "colnames" ].append( colname )
                    except KeyError:
                        t[ "colnames" ] = [ colname ]
                    sql[ c_id ] = colname
                    try:
                        h = "VARCHAR( {:d} )" . format( int( coltype ) )
                        h_def = " DEFAULT ''"
                    except ValueError:
                        if coltype.upper( ) == 'I':
                            h = "INT"
                            h_def = " DEFAULT 0"
                        elif coltype.upper( ) == 'T':
                            h = "TEXT"
                            h_def = " DEFAULT ''"
                        else:
                            h = "TEXT"
                            h_def = " DEFAULT ''"
                    # UNIQUE and INDEX keys are redundant, with UNIQUE taking priority
                    if "u" in colprops:  # add a UNIQUE KEY statement for this column
                        t[ "UNIQUE" ].append( c_id )
                    elif "i" in colprops:  # add an INDEX statement for this column
                        t[ "INDEX" ].append( c_id )
                    if "n" in colprops:  # add a NOT NULL statement for this colum
                        h += " NOT NULL" + h_def
                    #if i < len( t[ "depends" ] ):  # set up a matched temporary column if dependencies need to be resolved
                    #    # The temporary column gets the correct type assignment while the permanent one gets TEXT (so it will accomodate any value).
                    #    # Data is read into the permanent column in the beginning while the temporary column takes over after mapping.
                    #    sql[ c_id + _tmp ] = colname + _tmp
                    #    sql[ c_id + _tmp + _fmt ] = h
                    #    sql[ c_id + _fmt ] = "VARCHAR( 200 )"  # 'TEXT' does not work (easily) when keys are also set up for this column
                    #    t[ "col_ids_with_tmps" ].append( c_id + _tmp )
                    #else:
                    #    # no dependency -> no temporary column -> everything is fine
                    #    sql[ c_id + _fmt ] = h

                # decide whether this database participates in network building
                if "query" in t:
                    try:
                        t[ "query" ] = bool( int( t[ "query" ] ) )
                    except ValueError:
                        t[ "query" ] = t[ "query" ] in ( "yes", "y" )
                    if t[ "query" ]:
                        result[ "query_dbs" ].add( t[ "name" ] )
                #else:
                #    result[ "query_dbs" ].add( t[ "name" ] )

            except ( ConfigParser.NoOptionError, KeyError ):
                self._extalert( "Entries for config file section {!r} incomplete." . format( table ) )
        #print( dbt[ "interactions" ] ); sys.exit( )
        result[ "tables" ].update( {
            "taxa": {
                "update" : 1,
                "update-order" : 1,
                "name" : "Taxa",
                "tablename" : "Taxa",
                "sep" : "|",
                "write" : "sql",
                'sources': 'nodes.dmp.* names.dmp.*',
                'paths': ['./data/db_files/taxonomy/nodes.dmp.*', './data/db_files/taxonomy/names.dmp.*'],
                'dir': './data/db_files/taxonomy',
            },
            "taxon_aliases": {
                "update" : 1,
                "update-order" : 2,
                "name" : "Taxon_aliases",
                "tablename" : "Taxon_aliases",
                "sep" : "|",
                "write" : "sql",
                'sources': 'nodes.dmp.* names.dmp.*',
                'paths': ['./data/db_files/taxonomy/nodes.dmp.*', './data/db_files/taxonomy/names.dmp.*'],
                'dir': './data/db_files/taxonomy',
            },
            "actors": {
                "update" : 1,
                "update-order" : 2,
                "name" : "Actors",
                "tablename" : "Actors",
                "sep" : "\t",
                "write" : "sql",
                'sources': 'gene_info.* gene2accession.* uniprot_idmapping.* miRNA.dat miRNA_aliases.txt gene.tab factor.tab drugs.csv',
                'paths': ['./data/db_files/NCBI/gene_info.*', './data/db_files/NCBI/gene2accession.*', './data/db_files/NCBI/uniprot_idmapping.*', './data/db_files/NCBI/miRNA.dat', './data/db_files/NCBI/miRNA_aliases.txt', './data/db_files/NCBI/gene.tab', './data/db_files/NCBI/factor.tab', './data/db_files/NCBI/drugs.csv' ],
                'dir': './data/db_files/NCBI',
            },
            "actor_aliases": {
                "update" : 1,
                "update-order" : 3,
                "name" : "Actor_aliases",
                "tablename" : "Actor_aliases",
                "sep" : "\t",
                "write" : "sql",
                'sources': 'gene_info.* gene2accession.* uniprot_idmapping.* miRNA.dat miRNA_aliases.txt gene.tab factor.tab drugs.csv',
                'paths': ['./data/db_files/NCBI/gene_info.*', './data/db_files/NCBI/gene2accession.*', './data/db_files/NCBI/uniprot_idmapping.*', './data/db_files/NCBI/miRNA.dat', './data/db_files/NCBI/miRNA_aliases.txt', './data/db_files/NCBI/gene.tab', './data/db_files/NCBI/factor.tab', './data/db_files/NCBI/drugs.csv' ],
                'dir': './data/db_files/NCBI',
            },
            "actor_xrefs": {
                "update" : 1,
                "update-order" : 1,
                "name" : "Actor_xrefs",
                "tablename" : "Actor_xrefs",
                "sep" : "\t",
                "write" : "sql",
                'sources': 'MIRIAM.xml',
                'paths': ['./data/db_files/MIRIAM/MIRIAM.xml' ],
                'dir': './data/db_files/MIRIAM',
            },
            "actor_roles": {
                "update" : 1,
                "update-order" : 1,
                "name" : "Actor_roles",
                "tablename" : "Actor_roles",
                "sep" : "\t",
                "write" : "sql",
                'sources': 'gene_info.* gene2accession.* uniprot_idmapping.* miRNA.dat miRNA_aliases.txt gene.tab factor.tab drugs.csv',
                'paths': ['./data/db_files/NCBI/gene_info.*', './data/db_files/NCBI/gene2accession.*', './data/db_files/NCBI/uniprot_idmapping.*', './data/db_files/NCBI/miRNA.dat', './data/db_files/NCBI/miRNA_aliases.txt', './data/db_files/NCBI/gene.tab', './data/db_files/NCBI/factor.tab', './data/db_files/NCBI/drugs.csv' ],
                'dir': './data/db_files/NCBI',
            },
            "complexes": {
                "update" : 1,
                "update-order" : 3,
                "name" : "Complexes",
                "tablename" : "Complexes",
                "sep" : "\t",
                "write" : "sql",
                'sources': 'gene_info.* gene2accession.* uniprot_idmapping.* miRNA.dat miRNA_aliases.txt gene.tab factor.tab drugs.csv',
                'paths': ['./data/db_files/NCBI/gene_info.*', './data/db_files/NCBI/gene2accession.*', './data/db_files/NCBI/uniprot_idmapping.*', './data/db_files/NCBI/miRNA.dat', './data/db_files/NCBI/miRNA_aliases.txt', './data/db_files/NCBI/gene.tab', './data/db_files/NCBI/factor.tab', './data/db_files/NCBI/drugs.csv' ],
                'dir': './data/db_files/NCBI',
            },
        } )

        # parse dependencies per table
        # Dependencies are represented by nested lists: list of dependent columns -> list of relevant databases -> list of database couplings (cf. setup.cfg)
        #for t in result[ "tables" ].values( ):
        #    t[ "depends" ] = [ k.split( ":" ) for k in t[ "depends" ] ]  # separate database dependencies per column
        #    for i in xrange( len( t[ "depends" ] ) ):
        #        t[ "depends" ][ i ] = [ result[ "db_mapper" ][ l ] for l in t[ "depends" ][ i ] ]  # map identifier to table name

        #        # TODO: determine if database couplings will be implemented
        #        #t[ "depends" ][ i ] = [ k.split( ">>" ) for k in t[ "depends" ][ i ] ]  # separate database couplings per dependency
        #        #for j in xrange( len( t[ "depends" ][ i ] ) ):
        #        #    t[ "depends" ][ i ][ j ] = [ result[ "db_mapper" ][ l ] for l in t[ "depends" ][ i ][ j ] ]  # map identifier to table name

        # parse update order
        m = max( [ t[ "update" ] for t in result[ "tables" ].values( ) if type( t[ "update" ] ) == int ] ) + 1
        for t in result[ "tables" ].values( ):
            if t[ "update" ] in ( "yes" , "y" ):
                t[ "update" ] = m

        # graph settings
        nodes = ( "_fill", "_shape", "_label", "_height", "_width" )
        edges = ( "_color", "_source", "_target", "_width" )
        result[ "graphics" ] = dict( )
        moltypes = {  k.lower( ) for k in result[ "moltypes" ] | set( [ "seed" ] )  }
        for key in [ "xgmml" ]:
            result[ "graphics" ][ key ] = dict( conf.items( key ) )
            g = result[ "graphics" ][ key ]  # alias dictionary for easier access
            for k, v in g.items( ):
                if k == "extension":
                    continue
                a = v.split( )
                if len( a ) == 5:  # we assume a node entry
                    if k not in moltypes:
                        self._alert( "Option {!r} in section {!r} of the configuration file has no 'moltype' equivalent (cf. section 'miRNAworkbench')." . format( k, key ) )
                        continue
                    #kl = k.lower( )  # makes writing in the code easier NOTE: gotcha - option names are automatically lower-cased by ConfigParser?!?
                    try:
                        for i in ( 0, 2 ):  # fill color and label color
                            try:
                                g[ k + nodes[ i ] ] = "#{:06x}" . format( int( a[ i ], 16 ) )
                            except ValueError:
                                try:
                                    g[ k + nodes[ i ] ] = name_to_hex( a[ i ] )
                                except:
                                    self._extalert( "Value for option {!r} in section {!r} of the configuration file is invalid." . format( k, key ) )
                                    sys.exit( 1 )
                        for i in ( 3, 4 ):  # height and width
                            g[ k + nodes[ i ] ] = float( a[ i ] )
                        for i in ( 1, ):  # shape
                            g[ k + nodes[ i ] ] = a[ i ]
                    except:
                        self._extalert( "Value for option {!r} in section {!r} is invalid." . format( k, key ) )
                        sys.exit( 1 )
                elif len( a ) == 4:  # we assume an edge entry
                    try:
                        for i in ( 0, ):
                            g[ k + edges[ i ] ] = "#{:06x}" . format( int( a[ i ], 16 ) )
                        for i in ( 1, 2, 3 ):
                            g[ k + edges[ i ] ] = int( a[ i ] )
                    except:
                        self._extalert( "Value for option {!r} in section {!r} of the configuration file is invalid." . format( k, key ) )
                        sys.exit( 1 )
                else:  # unknown entry format
                    self._alert( "Unknown entry format for key {!r} in section {!r} of the configuration file." . format( k, key ) )

        return result


    def _add_to_conf( self, name, item ):
        self._conf_dict[ name ] = item


    def _get_table_list( self ):
        try:
            return sorted( self._conf_dict[ "known_dbs" ].values( ) )
        except KeyError:
            return [ ]


    def getConfAll( self ):
        return self._conf_dict


    def getConfItem( self, item ):
        """try to return the requested item"""
        try:
            return self._conf_dict[ item ]
        except KeyError:
            self._extalert( "No configuration item with name {!r}." . format( item ) )
            return None


    def getSQLMapping( self ):
        """return the mapping dictionary for SQL queries"""
        return self._conf_dict[ "sql_mapper" ]


    def getHardcodedTables( self ):
        """return the tables that are not up for user configuration"""
        return self._config.get( "database", "hardcoded" )


    def getTablesConf( self ):
        """return all table configuration"""
        return self._conf_dict[ "tables" ]


    def getTableConf( self, tname ):
        """return the dictionary of the given table"""
        try:
            return self._conf_dict[ "tables" ][ tname.lower( ) ]
        except:
            self._extalert( "No table with name {!r} in configuration." . format( tname ) )
            return None


    def getGraphicsConf( self, gformat ):
        """return the dictionary of the given table"""
        try:
            return self._conf_dict[ "graphics" ][ gformat ]
        except:
            self._extalert( "No table with name {!r} in configuration." . format( gformat ) )
            return None


    def _commit_db( self ):
        """induce the necessary SQL command to commit all executed INSERT, UPDATE, and DELETE statements"""
        self._serv.commit( )


