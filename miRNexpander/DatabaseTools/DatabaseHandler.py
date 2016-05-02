#!/usr/bin/env python

# DatabaseHandler.py
# Created by Martin Eberhardt ( martin.eberhardt@uk-erlangen.de )
# This file defines a class for working with the database:
# the initial setup, queries and updates my be performed.

from __future__ import print_function

import sys, os, re, tempfile
import math  # for log10 and ceil
import codecs  # for file objects ('open( )') with utf-8 writer support
from Bio import SeqIO
from operator import itemgetter
from glob import glob
from getpass import getpass
import MySQLdb
from xml.dom.minidom import parse
import time
from miRNexpander.DatabaseTools.DatabaseConnector import DatabaseConnector
from miRNexpander.DatabaseTools.SpeciesChecker import get_latin_abbr, get_scientific_name


class DatabaseHandler( DatabaseConnector ):
    """setup, update and query the database"""


    def __init__( self ):
        DatabaseConnector.__init__( self )

        self._setup_order = [ "DROP", "CREATE", "INSERT" ]  # the usual query order for database setups/updates
        self._db_input_methods = {
                "Actors" . lower( ) : self._read_genes_file,
                "Actor_aliases" . lower( ) : self._read_genes_file,
                "Actor_roles" . lower( ) : self._read_genes_file,
                "Actor_xrefs" . lower( ) : self._read_MIRIAM_file,
                "Complexes" . lower( ) : self._read_genes_file,
                "Drugbank" . lower( ) : self._read_Drugbank_file,
                "gene_xrefs" . lower( ) : self._read_genes_file,
                "HPRD" . lower( ) : self._read_HPRD_file,
                "HTRIdb" . lower( ) : self._read_HTRIdb_file,
                "miRecords" . lower( ) : self._read_miRecords_file,
                "miRTarBase" . lower( ) : self._read_miRTarBase_file,
                "miRTarBase6" . lower( ) : self._read_miRTarBase_file,
                "RegPhos" . lower( ) : self._read_RegPhos_file,
                "starBase" . lower( ) : self._read_starBase_file,
                "STRING" . lower( ) : self._read_STRING_file,
                "TarBase" . lower( ) : self._read_TarBase_file,
                "taxa" . lower( ) : self._unify_taxonomy,
                "taxon_aliases" . lower( ) : self._unify_taxonomy,
                "TRANSFAC" . lower( ) : self._read_TRANSFAC_file,
                "TransmiR" . lower( ) : self._read_TransmiR_file,
                }
        self._allowed_species = None
        self._db_restriction = None

        self.unknown_entity = int( 1e9 )  # located at border between genes and complexes; MySQL int is 4 bytes signed, so max is 2.147e9

        self._reg_type = {
                            "activation" : re.compile( "(?i)^activation|^induction|up|positive|double[- ]?negative" ),
                            "inhibition" : re.compile( "(?i)^inhibition|^repression|down|(?!double[- ]?)negative" ),
                        }
        self._key_replace = re.compile( "[^a-z0-9.]" )

        from MySQLdb.converters import conversions #as self._mysql_escape_conversions
        from _mysql import escape #as self._mysql_escape
        self._mysql_escape_conversions = conversions
        self._mysql_escape = escape
        del conversions, escape


    def connect( self, conf ):
        """add some method definitions to the configuration"""

        res = DatabaseConnector.connect( self, conf )
        self._db_restriction = self.getConfItem( "query_dbs" )
        self._add_to_conf( "db_input_methods", self._db_input_methods )
        return res


############################################################
#### SETUP/UPDATE METHODS                               ####
############################################################


    ############################################################
    #### PUBLIC                                             ####
    ############################################################

    def update_all( self ):
        """reread base tables and all interaction tables from database files"""

        # use fixed database schema and follow these steps:
        # - drop whole database by
        # read and load taxa, taxon_aliases, keep aliases for translation
        # read and load and link Actors, Actor_xrefs, Actor_aliases, Actor_Roles, Complexes, keep aliases and xrefs for translation
        # - populate gene.info first, then the rest (gene2accession, miRBase, factor.tab, etc)
        # read and load interaction tables by doing in-memory translation in python
        self._conf_dict[ "known_dbs" ] = {
                "taxa": "Taxa",
                "taxon_aliases" : "Taxon_aliases",
                "actors" : "Actors",
                "actor_xrefs" : "Actor_xrefs",
                "actor_aliases" : "Actor_aliases",
                "actor_roles" : "Actor_roles",
                "complexes" : "Complexes",
                "mirbase" : "miRBase",
        }
        # empty the database by redoing the setup; this drops the whole database schema and is very fast
        user = self._config.get( 'database', 'user' )
        db = self._config.get( 'database', 'database' )
        prompt = """You have requested to update the whole database. This may take minutes to several hours to complete.
You will have to wait for the update to finish before you can run queries
WARNING: Before anything else, all database content will be dropped if you enter the correct password!
To proceed, enter the database password for MySQL user {}: """ . format( user )
        if os.system( "mysql -u {} --password={} {} < database_schema.sql" . format( user, getpass( prompt ), db ) ):
            self._spill( "There was a problem while emptying the database. Please do a manual setup using database_schema.sql!" )
            sys.exit( )
        self._spill( "Database was emptied - starting update process..." )
        for tables in [ "taxa", "taxon_aliases", "actor_xrefs" ], [ "actor_roles", "actors", "actor_aliases", "complexes", "mirbase" ]:
            self.update_db( tables )


    def update_db( self, table_list, target_list = [ ], source_list = [ ], always_update = False ):
        """reload the specified databases, write SQL to the specified files, use the specified input files"""

        # 1st condition line checks: update the listed databases, or all of them (if True)
        # 2nd condition line checks: write to the listed sql files, or the standard ones (if True)
        # 3rd condition line checks: read from the listed data files, or use standard ones (if None)
        if type( table_list ) != list and table_list != True \
            or type( target_list ) != list and target_list != True  \
            or type( source_list ) != list and source_list != None:
            self._extalert( "This method takes at most three lists as arguments." )
            return 1

        if table_list == True:
            tables = set( [  t for t in self._get_table_list( ) if self.getTableConf( t )[ "update" ] not in ( 0, "n", "no" )  ] )
            unspecific = True
        else:
            tables = set( table_list )
            unspecific = False

        # safety valve
        #tables -= set( [ "actor_xrefs", "actors", "actor_aliases", "taxa", "taxon_aliases" ] )

        # sort tables according to their configured update order
        table_list = dict( [  ( t, self.getTableConf( t )[ "update-order" ] ) for t in tables  ] )
        table_list = zip( *sorted( table_list.iteritems( ), key = itemgetter( 1 ) ) )[ 0 ]
        lowered_set = { t.lower( ) for t in table_list }

        #print( "List: ", table_list )
        for i in xrange( len( table_list ) ):
            table = table_list[ i ]
            tl = table.lower( )
            if tl not in lowered_set:  # prevents repeated processing (s.b.)
                continue

            #print( "Set: ", lowered_set )
            t_conf = self.getTableConf( tl )
            try:
                table_id = t_conf[ "sources" ]
            except ( TypeError, KeyError ):
                self._alert( "Warning: Database {!r} is unknown and will be ignored." . format( table ) )
                continue

            table_name = t_conf[ "name" ]

            if not always_update:
                build = t_conf[ "update" ]
                if unspecific and not build:
                    self._spill( "Updates for database {!r} disabled in configuration, skipping." . format( table ) )
                    continue

            self._spill( "Processing input for {}..." . format( table_name ) )
            try:
                db_files = source_list[ i ].split( "," )  # fails if source_list[ i ] is non-existing or None
            except ( IndexError, TypeError ):
                db_files = t_conf[ "paths" ]

            sep = t_conf[ "sep" ]
            if sep == r"\t":
                sep = "\t"

            try:
                input_method = self.getConfItem( "db_input_methods" )[ tl ]
            except KeyError:
                self._alert( "Warning: No method definition for database {!r}, skipping update." . format( table ) )
                continue

            inserts = input_method( db_files, sep )

            try:
                write_file = target_list[ i ]
            except ( TypeError, IndexError ):
                if target_list == True:
                    write_file = inserts.keys( )
                else:
                    write_file = False

            ### internal subroutine ###
            def _insert_template( out_type, elements, row_length = None ):
                """return a structured string of MySQL-escaped values"""

                count = len( elements )
                if count < 1:
                    return [ "" ]

                if type( row_length ) is int:
                    # separate elements into multiple rows
                    row_length = min( row_length, count )  # make sure to insert at most count spaceholders
                    rows = count // row_length
                else:
                    row_length = count
                    rows = 1

                if out_type == "sql":
                    s = [ ]
                    row_limit_per_query = int( 8e4 )
                    c = 0
                    while c < rows:
                        rem = min( row_limit_per_query, rows - c )
                        h = ",\n\t" . join( [ "( " + ", " . join( [ "{}" ] * row_length ) + " )" ] * rem ) + ";\n\n"
                        s.append( h.format( *self._mysql_escape( elements[ c * row_length : ( c + rem ) * row_length ], self._mysql_escape_conversions ) ) )
                        c += row_limit_per_query
                elif out_type == "txt":
                    s = [ "\n" . join( [ "\t" . join( [ "{}" ] * row_length ) ] * rows ) + "\n" ]
                    s[ 0 ] = s[ 0 ].format( *self._mysql_escape( elements, self._mysql_escape_conversions ) )
                else:
                    s = [ "" ]
                return s

            ### end of internal subroutine ###

            ( queries, elements_per_row ) = self._prepare_setup_frame( inserts.keys( ) )
            #print( table_list, inserts.keys( ), queries.keys( ), elements_per_row.keys( ) )
            for table in table_list:
                if table.lower( ) not in {  q.lower( ) for q in queries  }:  # prevents another table that is read from the same file from updating
                    continue
                lowered_set.remove( table.lower( ) )  # prevents repeated processing (s.a.)

                # convolution revisited: turn strings into an one-element iterable of the corresponding string
                # -> This all became necessary because MySQL's query_cache_size parameter limits the length of query strings.
                # -> Instead of INSERTing 1 000 000 rows in one query, break up the inserts and run multiple queries.
                # -> Since later scripts expect a unitary data structure, all queries are turned into iterables here.
                for i in ( 0, 1 ):
                    queries[ table ][ self._setup_order[ i ] ] = ( queries[ table ][ self._setup_order[ i ] ], )
                # add inserts
                ins_stmt = queries[ table ][ self._setup_order[ 2 ] ]
                queries[ table ][ self._setup_order[ 2 ] ] = [ ]
                #for ins in _insert_template( t_conf[ "write" ], inserts[ table ], elements_per_row[ table ] ):
                for ins in _insert_template( "sql", inserts[ table ], elements_per_row[ table ] ):
                    queries[ table ][ self._setup_order[ 2 ] ].append( ins_stmt + ins )
#            tmp = tempfile.NamedTemporaryFile( )
#            tmp.write( queries[ -1 ] )
#            LOAD DATA LOCAL INFILE "/tmp/taxa.sql" REPLACE INTO TABLE taxa IGNORE 1 LINES
#            Need to start mysql client with --local-infile=1 for the above to work.

                self._spill( "Processing table setup for {}..." . format( table ) )
                with codecs.open( "data/sql/" + table + ".sql", "w", "utf-8" ) as f:
                    for l in queries[ table ][ "INSERT" ]:
                        f.write( l )
                #self._execute_setup_queries( t_conf[ "write" ], queries[ table ] )
                if table != "actor_aliases":
                    self._execute_setup_queries( "sql", queries[ table ] )


            #if table not in self.getHardcodedTables( ):
            #    self._spill( "\t->Resolving IDs..." )
            #    self._resolve_table_dependencies( queries )
            self._spill( "...Table setup finished." )
            self._store_setup_data( t_conf[ "write" ], queries, write_file )


    ############################################################
    #### GENERAL AND NOT SO PUBLIC                          ####
    ############################################################

    def _prepare_setup_frame( self, db_ids, out_type = "sql", required_additional_tables = [ ], out = None ):
        """build the fundamental creation query for the specified table"""

        ## We use the string method 'format' together with the above dictionary to insert
        ## values from the configuration file. The advantages are cross-module consistency
        ## and easy adaptation through edits to the configuration file.
        q = self.getSQLMapping( )

        if type( db_ids ) not in ( tuple, list ):
            db_ids = [ db_ids ]

        if set( [ "taxa", "taxon_aliases" ] ) & set( db_ids ):
            queries = {
                "taxa" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Taxa VALUES\n\t"
                },
                "taxon_aliases" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Taxon_aliases VALUES\n\t"
                },
            }
            elements = {
                "taxa" : 3,
                "taxon_aliases" : 2,
            }
            return( queries, elements )
        elif set( [ "actors", "actor_roles", "actor_xrefs", "actor_aliases", "complexes", "mirbase" ] ) & set( db_ids ):
            queries = {
                "actor_roles" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Actor_roles VALUES\n\t",
                },
                "actors" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Actors VALUES\n\t",
                },
                "actor_aliases" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Actor_aliases VALUES\n\t",
                },
                "actor_xrefs" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Actor_xrefs VALUES\n\t",
                },
                "complexes" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO Complexes VALUES\n\t",
                },
                "mirbase" : {
                    "CREATE" : "SHOW TABLES;",
                    "DROP" : "SHOW TABLES;",
                    "INSERT" : "INSERT INTO miRBase VALUES\n\t",
                },
            }
            elements = {
                "actors" : 4,
                "actor_aliases" : 3,
                "actor_xrefs" : 7,
                "actor_roles" : 2,
                "complexes" : 4,
                "mirbase" : 5,
            }
            return( queries, elements )

        queries, elements = { }, { }
        for table_id in db_ids:
            if type( table_id ) != str: continue

            # extract information for this table
            if table_id in self.getHardcodedTables( ):
                t = dict( name = table_id, columns = self._sql( "SHOW COLUMNS IN {};" . format( table_id ) ) )
                delete_table_stmt = "DELETE FROM {};\n" . format( t[ "name" ] )
                create_table_stmt = "SHOW TABLES;\n"  # stand in since MySQL does not accept empty commands
                hardcoded = True
            else:
                t = self.getTableConf( table_id )
                delete_table_stmt = "DELETE FROM {};\n" . format( t[ "name" ] )
                create_table_stmt = "SHOW TABLES;\n"  # stand in since MySQL does not accept empty commands

            if out_type == "sql":
                insert_scaffold = "INSERT INTO {} VALUES\n\t" . format( t[ "name" ] )
            elif out_type == "txt":
                insert_scaffold = ""
            else:
                self._alert( "Unknown file type for output of table {!r}: {!r}." . format( table_id, out_type ) )
                continue

            queries[ table_id ] = {
                                self._setup_order[ 0 ] : delete_table_stmt,
                                self._setup_order[ 1 ] : create_table_stmt,
                                self._setup_order[ 2 ] : insert_scaffold,
                                }
            elements[ table_id ] = 5  # NOTE: needs to be read from the table configuration !!!

        return ( queries, elements )


    def _execute_setup_queries( self, out_type, queries ):
        """execute the queries and commit the database"""

        if out_type == "sql":
            automatic = [  self._setup_order[ i ] for i in ( 0, 2 )  ]
        elif out_type == "txt":
            automatic = self._setup_order[ : -1 ]

        for query_type in automatic:
            self._spill( "\t-> Executing " + query_type + "..." )
            for q in queries[ query_type ]:
                try:
                    self._cursor.execute( q )
                except MySQLdb.Warning:
                    pass
                except MySQLdb.Error as e:
                    self._alert( "SQL QUERY Error.\nDescription: {}\nQuery started with:\n{!r}\n" . format( e, q[ 0:76 ] ) )
                    continue

        if out_type == "txt":
            self._spill( "\t-> Executing " + self._setup_order[ -1 ] + "..." )

        self._commit_db( )  # this is immensely important to call, otherwise all SQL INSERT, UPDATE, and DELETE statements will not be finalized


    def _store_setup_data( self, out_type, queries, targets ):
        """write the given queries to a permanent buffer (file)"""

        if type( targets ) == str:
            target = targets
            try:
                sql_file = open( target, "w" )
                sql_file.write( "\n" . join( [ "" . join( queries[ table ][ qtype ] ) for table in queries for qtype in self._setup_order ] ) )
                sql_file.close( )
            except IOError as e:
                self._extalert( e.strerror )
        elif type( targets ) in ( dict, list, set ):
            for t in targets:
                try:
                    target = os.path.join( self.getConfItem( "sql_path" ), self.getTableConf( t )[ "name" ] + "." + out_type )
                    sql_file = open( target, "w" )
                    sql_file.write( "\n" . join( [ "" . join( queries[ t ][ q ] ) for q in self._setup_order ] ) )
                    sql_file.close( )
                except IOError as e:
                    self._alert( e.strerror )
                except ( KeyError, AttributeError ):
                    self._extalert( "Mismatching keys between 'queries' and 'targets'.\n\tqueries: ", queries.keys( ), "\n\ttargets: ", targets.keys( )  )
        elif targets == True:  # temporary file management not implemented yet
            pass


    def _build_std_sql_target_list( self ):
        """generates a list of files to write SQL statements to"""

        path = self.getConfItem( "sql_path" )
        result = [ ]
        for db in self._get_table_list( ):
            result.append( os.path.join( path, "{}.sql" . format( db ) ) )
        return result


    ############################################################
    #### MIRIAM annotation database                         ####
    ############################################################

    ### encapsulated method for reading data from an XML file ###
    def _read_MIRIAM_file( self, db_files, sep = None ):
        """read identifier information from MIRIAM file"""

        if type( db_files ) not in ( list, set, dict ):
            self._extalert( "Wrong type for argument 1: expected iterable." )
            return { }

        MIRIAM = dict( )  # result variable
        for db_file in db_files:
            try:
                xml_object = parse( db_file )
            except IOError:
                self._spill( "Unable to read input file {!r}, aborting." . format( db_file ) )
                continue
            except xml.parsers.expat.ExpatError as e:
                self._spill( "There was a parsing error, please check if {!r} is a valid XML document!" . format( db_file ) )
                self._spill( "Error message:\n\t" + e.str )
                continue

            columns = set( )
            for x in xml_object.getElementsByTagName( "datatype" ):
                d = dict( [ ( a, x.getAttribute( a ) ) for a in x.attributes.keys( ) ] )
                d.update( [ ( k, x.getElementsByTagName( k )[ 0 ].childNodes[ 0 ].data ) for k in ( "name", "namespace", "definition" ) ] )
                for u in x.getElementsByTagName( "uri" ):
                    if not u.getAttribute( "deprecated" ):
                        d[ u.getAttribute( "type" ) ] = u.childNodes[ 0 ].data

                MIRIAM[ d[ "id" ] ] = d
                columns |= set( d )

        insert_values = list( )
        i = 0
        for k, v in sorted( MIRIAM.items( ) ):
            i += 1
            k = k.encode( "ascii", "xmlcharrefreplace" )
            # utf-8 replacement is (hopefully) temporary
            for x, y in v.items( ):
                v[ x ] = y.encode( "ascii", "xmlcharrefreplace" )
            v.update( [ ( col, "" ) for col in columns - set( v ) ] )  # generate empty values for potentially missing attributes
            insert_values.extend( [ i, k, v[ "name" ], v[ "namespace" ], v[ "definition" ], v[ "URN" ], v[ "URL" ] ] )

        return {  "actor_xrefs" : insert_values  }


    ############################################################
    #### NCBI Taxonomy Nodes                                ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_NCBI_taxnode_files( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        if type( db_files ) not in ( list, set, dict ):
            self._extalert( "Wrong type for argument 1: expected iterable." )
            return { }

        taxa = dict( )  # result variable
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._alert( e.strerror + ": {}" . format ( db_file ) )
                return taxa

            for record in stream:
                line = record.split( sep )
                try:
                    tax_id = int( line[ 0 ].strip( ) )
                    parent_id = int( line[ 1 ].strip( ) )
                    taxa[ tax_id ] = {  "parent" : parent_id  }

                except ( IndexError, ValueError ):
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        return taxa


    ############################################################
    #### NCBI Taxonomy Names                                ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_NCBI_taxname_files( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        taxa = dict( )  # result variable
        if type( db_files ) not in ( list, set, dict ):
            self._extalert( "Wrong type for argument 1: expected iterable." )
            return taxa

        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._alert( e.strerror + ": {}" . format ( db_file ) )
                return taxa

            for record in stream:
                line = record.split( sep )
                try:
                    tax_id = int( line[ 0 ].strip( ) )
                    alias = line[ 2 ].strip( )
                    if alias == "":  # most unique names will not be set
                        alias = line[ 1 ].strip( )
                    alias_type = line[ 3 ].strip( )
                    if alias_type in ( "scientific name", ):#"common name" ):
                        try:
                            taxa[ tax_id ][ alias_type ] = alias
                        except KeyError:
                            taxa[ tax_id ] = {  alias_type : alias  }
                    else:
                        try:
                            taxa[ tax_id ][ "aliases" ] |= set( [ alias ] )
                        except KeyError:
                            try:
                                taxa[ tax_id ][ "aliases" ] = set( [ alias ] )
                            except KeyError:
                                taxa[ tax_id ] = {  "aliases" : set( [ alias ] )  }

                except ( IndexError, ValueError ):
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        return taxa



    ############################################################
    #### NCBI Taxonomy Unifier                              ####
    ############################################################

    def _unify_taxonomy( self, dummy1, dummy2 ):
        """merge the data structures returned by two other methods"""

        taxo = self.getTableConf( "taxa" )

        nodes_files, names_files = map( glob, taxo[ "paths" ] )
        names_sep = nodes_sep = "|"

        d1 = self._read_NCBI_taxnode_files( nodes_files, nodes_sep )
        d2 = self._read_NCBI_taxname_files( names_files, names_sep )

        insert_values = {  "taxa" : [ ], "taxon_aliases" : [ ]  }

        for k in d1:
            try:
                node = d1[ k ]
                names = d2[ k ]
            except KeyError:  # means there was no entry in d2 for k, so no names were found
                continue

            insert_values[ "taxon_aliases" ].extend( [ k, str( k ), k, names[ "scientific name" ] ] )
            insert_values[ "taxa" ].extend( [ node[ "parent" ], k, names[ "scientific name" ] ] )
            try:
                for a in names[ "aliases" ]:
                    insert_values[ "taxon_aliases" ].extend( [ k, a ] )
            except KeyError:
                pass

        return insert_values


    ############################################################
    #### genes                                              ####
    ############################################################

    ### encapsulated method for compiling data from database files ###
    def _read_genes_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        monicker = "actor"
        tAnnot = monicker + "s"
        tAlias = monicker + "_aliases"
        tXrefs = monicker + "_xrefs"
        tRoles = monicker + "_roles"
        tCompl = "complexes"
        insert_values = {
                tAnnot : [ ],
                tAlias : [ ],
                tCompl : [ ],
        }

        miriam = dict( self._sql( "SELECT `namespace`, `x_id` FROM `Actor_xrefs`" ) )
        # calculate lowest index for additional databases to be inserted (MIRIAM does not account for TRANSFAC, Vega, etc)
        lmiriam = len( miriam )
        nm = "NEW:{:08d}"
        bmiriam = 10 ** int( math.ceil( math.log10( lmiriam ) ) )
        nmiriam = bmiriam + 1
        miriam[ "synonym" ] = nmiriam
        new_miriams = { ( nmiriam, nm.format( nmiriam ), "Synonyms", "synonym", "General synonyms, often human-readable", "urn:new:synonym", "" ) }
        nmiriam += 1

        # add aliases for some namespaces that are designated differently in the files to be read
        miriam.update( {
            "MIM" . lower( ) : miriam[ "omim" ],
        } )

        no_gene = self.unknown_entity

        annot = { }
        complexes = { }
        drugs = { }
        xrefs = { no_gene : set( [ ( miriam[ "ncbigene" ], "0" ) ] ) }
        entrez_mapper = { }
        aliased_mirnas = { }
        symbols = { }  # necessary to translate between gene_info and factor.tab
        uniprot = { }  # necessary to translate between gene_info and drugs.csv
        transfac = { }  # necessary to translate between gene.tab and factor.tab + site.tab
        sep = "\t"
        rem = set( [ "", "-" ] )
        roles = set( [ "protein", "miRNA", "gene", "transcript" ] )
        dbf = [ ]
        for f in map( glob, db_files ):
            dbf.extend( f )
        db_files = dbf
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._alert( "Failed to read {}" . format( db_file ) )
                continue
            self._spill( "\t-> Reading {}..." . format( os.path.abspath( db_file ) ) )
            fn = os.path.basename( db_file )

            if fn[ :9 ] == "gene_info":
                # add every single obscure and un-MIRIAMed database that NCBI loves to reference
                miriam[ "imgt/gene-db" ] = nmiriam
                new_miriams.add( ( nmiriam, nm.format( nmiriam ), "IMG/Gene-DB", "imgt.gene_db", "from NCBI's gene_info", "urn:new:imgt.gene_db", "http://www.imgt.org/IMGT_GENE-DB/GENElect?species=Homo+sapiens&query=2+" ) )
                nmiriam += 1

                miriam[ "mgi" ] = nmiriam
                new_miriams.add( ( nmiriam, nm.format( nmiriam ), "MGI", "mgi", "JAXLab's MGI mouse database", "urn:new:mgi", "http://www.informatics.jax.org/marker/" ) )
                nmiriam += 1

                miriam[ "vega" ] = nmiriam
                new_miriams.add( ( nmiriam, nm.format( nmiriam ), "VEGA", "vega", "from NCBI's gene_info", "urn:new:vega", "http://vega.sanger.ac.uk/id/" ) )
                nmiriam += 1

                index = len( annot )
                for record in stream:
                    line = record.split( sep )
                    try:
                        tax_id = line[ 0 ].strip( )
                        gene_id = int( line[ 1 ].strip( ) )
                        symbol = line[ 2 ].strip( )
                        desc = line[ 8 ].strip( )
                        info = {  "tax_id" : int( tax_id ), "symbol" : symbol, "desc" : desc   }#, "Entrez" : gene_id  }

                        synonyms = set( line[ 4 ].strip( ).split( "|" ) ) - rem
                        xref = set( zip( [ miriam[ "synonym" ] ] * len( synonyms ), synonyms ) )
                        xref.add( ( miriam[ "hgnc.symbol" ], symbol ), )
                        xref.add( ( miriam[ "ncbigene" ], str( gene_id ) ), )
                        symbols[ tax_id + ":" + symbol ] = gene_id
                        for dbxref in set( line[ 5 ].strip( ).split( "|" ) ) - rem:
                            entry = dbxref.split( ":" )
                            entry[ 0 ] = miriam[ entry[ 0 ].lower( ) ]
                            entry[ 1 ] = ":" . join( entry[ 1: ] )  # necessary to reconstitute identifiers that contain a colon, e.g. MGI markers
                            xref.add( tuple( entry ) )

                        annot[ gene_id ] = info
                        xrefs[ gene_id ] = xref

                    except ( IndexError, ValueError ):
                        self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )
                stream.close( )
            elif fn[ :8 ] == "gene.tab":
                #AC ID  SD  OS  DR  BS  RX
                try:
                    miriam[ "transfac" ]
                except KeyError:
                    miriam[ "transfac" ] = nmiriam
                    new_miriams.add( ( nmiriam, nm.format( nmiriam ), "TRANSFAC", "transfac", "TRANSFAC identifiers", "urn:new:transfac", "" ) )
                    nmiriam += 1
                alias_pos = ( 0, miriam[ "transfac" ] ),
                for record in stream:
                    line = record.split( "\t" )
                    if line[ 3 ] not in ( "Homo sapiens", "Mus musculus" ):
                        continue
                    gene_id = None
                    aliases = [ ]
                    for r in line[ 4 ].split( "," ):
                        try:
                            gene_id = int( r )
                        except ValueError:
                            aliases.append( ( miriam[ "ensembl" ], r ) )
                    if not gene_id:
                        continue
                    aliases.append( ( miriam[ "transfac" ], line[ 0 ] ) )
                    transfac[ line[ 0 ] ] = gene_id
                    try:
                        xrefs[ gene_id ] |= set( aliases )
                    except KeyError:
                        pass
                stream.close( )
            elif fn[ :10 ] == "factor.tab":
                #AC  ID  FA  GE  OS  DR  BS  RX
                # empty GE means complex or derived factor (e.g. phosphorylation); translating those must rely on the factor symbol in column FA
                try:
                    miriam[ "transfac" ]
                except KeyError:
                    miriam[ "transfac" ] = nmiriam
                    new_miriams.add( ( nmiriam, nm.format( nmiriam ), "TRANSFAC", "transfac", "TRANSFAC identifiers", "urn:new:transfac", "" ) )
                    nmiriam += 1
                alias_pos = ( 0, miriam[ "transfac" ] ), ( 3, miriam[ "transfac" ] )
                unknown = set( )
                # define observed symbol mutilations using regular expressions
                removexps = [
                    # first order
                    re.compile( "^\(([^)]+)\)(\d+)$" ),  # e.g. (STAT3)2

                    # second order
                    re.compile( "^(.+)-isoform\d+[A-Za-z]*$" ),  # e.g. NF-kappaB-isoform1
                    re.compile( "^([^{]+)(\{[^}]+\})+$" ),  # e.g. STAT1{pS478}
                    re.compile( "^(.+)\([lI]\)$" ),  # e.g. Brn-3a(l)  (T24145)
                    re.compile( "^(.+) \(\d+\)$" ),  # e.g. AP-3 (1)
                    re.compile( "^(.+)\(-like\)$" ),  # e.g. NF-kappaB(-like)
                    re.compile( "^(.+)-(LIP|LAP|FL)$" ),  # e.g. C/EBPbeta-LAP
                    re.compile( "^(.+)-?([Ll]ong|[Ss]hort)$" ),  # e.g. GATA6short:GATA6long
                ]
                try:
                    tr_factors = dict( [ ( l[ 0 ], l[ 1 ] ) for l in [ l.lower( ).split( "\t" ) for l in open( os.path.dirname( db_file ) + "/TRANSFAC.unmappables" ) ] ] )
                except IOError:
                    self._spill( "No file found for mapping dirty TRANSFAC entries: {}." . format( os.path.dirname( db_file ) ) )
                    tr_factors = dict( )

                i = no_gene
                for record in stream:
                    line = record.split( "\t" )
                    if not line[ 7 ] or line[ 4 ] not in ( "Homo sapiens", "Mus musculus", "Mammalia" ):
                        continue
                    if line[ 4  ][ 0 ] == "H":
                        tax_prefix = "9606:"
                        case_convert = str.upper
                    else:
                        tax_prefix = "10090:"
                        case_convert = str.capitalize
                    # try to match the TRANSFAC id (GXXXXXXX) with an existing entry from gene.tab
                    line[ 3 ] = line[ 3 ].strip( )
                    try:
                        xrefs[ transfac[ line[ 3 ] ] ].add( ( miriam[ "transfac" ], line[ 0 ].strip( ) ), )
                    except KeyError:  # emtpy GE
                        members = line[ 2 ].strip( ).split( ":" )
                        mem = [ ]
                        # column FA uses synonyms instead of official identifiers, this makes a lot of black magic necessary
                        ## step 1: resolve pesky "(ABC)2" strings by transforming them into nice and countable "ABC:ABC"
                        for m in members:
                            x = removexps[ 0 ].match( m )
                            try:
                                mem.extend( [ x.group( 1 ) ] * int( x.group( 2 ) ) )
                            except AttributeError:
                                mem.append( m )
                        members = mem[ : ]
                        for j in xrange( len( members ) ):
                            m = members[ j ]
                            # try to use the imported rules for translating the factors' names
                            try:
                                m = tr_factors[ m.lower( ) ]
                            except KeyError:
                                for r in removexps[ 1: ]:
                                    m = r.sub( r'\1', m )
                                m = m.lower( )
                                # replace latin particles with their first letter
                                for p in "alpha", "beta", "gamma", "delta", "epsilon", "kappa":
                                    m = m.replace( p, p[ 0 ] )
                                # remove hyphens and slashes, correct capitalization according to species entry
                                m = case_convert( m.replace( "-", "" ).replace( "/", "" ) )
                            else:
                                # capitalization according to species
                                if m == "":
                                    m = "unmappable_identifier"
                                m = case_convert( m )
                            members[ j ] = m
                        if len( members ) == 1:  # monomeric -> add factor ID to corresponding gene aliases
                            try:  # try to match the symbol with an existing entry from gene.info
                                xrefs[ symbols[ tax_prefix + members[ 0 ].strip( ) ] ].add( ( miriam[ "transfac" ], line[ 0 ].strip( ) ), )
                            except KeyError:  # mark alias as unknown
                                xrefs[ no_gene ].add( ( miriam[ "transfac" ], line[ 0 ].strip( ) ), )
                                unknown.add( "\t" . join( [ tax_prefix, mem[ 0 ].strip( ), members[ 0 ].strip( ), line[ 7 ].strip( ), line[ 0 ] ] ) )
                        else:  # multimeric -> create a complex entry and add the factor ID as alias
                            members = sorted( members )  # keep subunit order fixed to ease analyses like occurrence counting etc.
                            member_genes = set( )
                            symbol = ":" . join( members )
                            i += 1
                            complexes[ i ] = { "tax_id" : int( tax_prefix[ :-1 ] ), "symbol" : symbol, "desc" : "from TRANSFAC", "members" : dict( ) }
                            for j in xrange( len( members ) ):
                                m = members[ j ]
                                try:
                                    gene_id = symbols[ tax_prefix + m ]
                                    member_genes.add( gene_id )
                                except KeyError:
                                    gene_id = no_gene
                                    member_genes.add( 0 )
                                    unknown.add( "\t" . join( [ tax_prefix, mem[ j ].strip( ), m, line[ 7 ].strip( ), line[ 0 ] ] ) )
                                # count how often the member features in the complex
                                try:
                                    complexes[ i ][ "members" ][ ( gene_id, "protein" ) ] += 1
                                except KeyError:
                                    complexes[ i ][ "members" ][ ( gene_id, "protein" ) ] = 1
                            # insert a derived Entrez ID alias for the complex
                            xrefs[ i ] = set( [
                                ( miriam[ "transfac" ], line[ 0 ].strip( ) ),
                                ( miriam[ "synonym" ], line[ 2 ] ),
                                ( miriam[ "hgnc.symbol" ], symbol ),
                                ( miriam[ "ncbigene" ], "," . join( map( str, member_genes ) ) ),
                            ] )
                #[  print( l ) for l in unknown  ]
            elif fn[ :14 ] == "gene2accession":
                # the format is as follows: (first column has index 0):
                # column 3 and 4 refer to RNA accessions, with 3 being either EMBL or RefSeq
                # column 5 and 6 refer to protein accessions, with 5 being either Uniprot, EMBL, RefSeq, or RPF (last one has very few cases)
                # columns 13 and 14 refer to mature peptide accessions, but only about 160 are present in a file with 50 million entries
                alias_pos = ( 3, miriam[ "ensembl" ] ), ( 4, miriam[ "ncbigi" ] ), ( 5, miriam[ "uniprot" ] ), ( 6, miriam[ "ncbigi" ] ), ( 13, miriam[ "refseq" ] ), ( 14, miriam[ "ncbigi" ] )
                for record in stream:
                    line = record.split( sep )
                    try:
                        gene_id = int( line[ 1 ].strip( ) )
                    except (IndexError, ValueError ):
                        continue
                    try:
                        for p, ref_type in alias_pos:
                            if line[ p ] not in rem:
                                line[ p ] = line[ p ].strip( ).split( "." )[ 0 ]
                                if p == 5:
                                    if "_" in line[ p ]:
                                        ref_type = miriam[ "refseq" ]
                                    elif line[ p ][ 0 ] == "0":
                                        ref_type = miriam[ "ncbigi" ]
                                    elif len( line[ p ] ) > 6:
                                        ref_type = miriam[ "ensembl" ]
                                    else:
                                        uniprot[ line[ p ] ] = gene_id
                                elif p == 3 and "_" in line[ p ]:
                                    ref_type = miriam[ "refseq" ]
                                xrefs[ gene_id ].add( ( ref_type, line[ p ] ), )
                    except IndexError:
                        self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )
                stream.close( )
            elif fn[ :17 ] == "uniprot_idmapping":
                alias_pos = ( 0, miriam[ "uniprot" ] ), ( 3, miriam[ "refseq" ] ), ( 13, miriam[ "omim" ] ), ( 16, miriam[ "ensembl" ] ), ( 18, miriam[ "ensembl" ] ), ( 19, miriam[ "ensembl" ] ), ( 20, miriam[ "ensembl" ] )
                for record in stream:
                    line = record.split( sep )
                    try:
                        gene_id = int( line[ 2 ].strip( ) )
                    except ( IndexError, ValueError ):
                        continue
                    try:
                        uniprot[ line[ 0 ] ] = gene_id
                        for p, ref_type in alias_pos:
                            helper = {  e.split( "." )[ 0 ] for e in line[ p ].strip( ).split( "; " )  } - rem
                            xrefs[ gene_id ] |= set( zip( [ ref_type ] * len( helper ), helper ) )
                    except KeyError:
                        pass  # discontinued records
                    except ( ValueError, IndexError ):
                        self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )
                stream.close( )
            elif fn[ -9: ] == "miRNA.dat":
                stream.close( )
                mapped_mirnas = dict( )
                mirna_refs = dict( )
                # access all the data inside the .dat file
                for seq_record in SeqIO.parse( db_file, "embl" ):

                    if seq_record.name[ :3 ] not in ( "hsa", "mmu" ):  # NOTE: only include mouse and human for now
                        continue
                    d = dict( [ e.split( ":" )[ :2 ] for e in seq_record.dbxrefs ] )
                    try:
                        gene_id = int( d[ "ENTREZGENE" ] )
                    except ( KeyError, ValueError ):
                        continue

                    mapped_mirnas[ seq_record.id ] = gene_id
                    mirna_refs[ seq_record.id ] = dict( derived = set( ), PMIDs = set( ) )
                    for feature in seq_record.features:
                        # the acc number
                        try:
                            acc = feature.qualifiers[ 'accession' ][ 0 ]
                        except KeyError:  # means that this is another feature, like 'modified base'
                            continue
                        else:
                            mapped_mirnas[ acc ] = gene_id  # link between precursor miRNA and the corresponding gene
                            mirna_refs[ seq_record.id ][ "derived" ].add( acc )  # link between precursor and mature miRNA

                    for ref in seq_record.annotations[ "references" ]:
                        try:
                            mirna_refs[ seq_record.id ][ "PMIDs" ].add( ref.pubmed_id )
                        except AttributeError:
                            continue

            elif fn[ -17: ] == "miRNA_aliases.txt":
                border = re.compile( "\s*;\s*" )
                # read in aliases
                for line in stream:
                    # line format is "miRBase_Accession    alias1;alias2;...", ending on a semicolon
                    acc, aliases = line.split( )
                    if aliases[ :3 ] not in ( "mmu", "hsa" ):  # NOTE: only include mouse and human for now
                        continue
                    acc = acc.strip( )
                    iv = border.split( aliases )[ :-1 ]  # get names, discarding final value which will always be ''
                    iv.reverse( )  # puts the current name up front
                    iv.extend( [  n[ 4: ] for n in iv if n[ :4 ] in ( 'hsa-', 'mmu-' )  ] )  # map missing species prefixes to human and mouse per default

                    aliased_mirnas[ acc ] = set( zip( [ miriam[ "synonym" ] ] * len( iv ), iv ) )
                stream.close( )
            elif fn[ -9: ] == "drugs.csv":  # Drugbank
                i = 2 * no_gene  # start counting from 2e9
                alias_pos = ( ( 0, miriam[ "drugbank" ] ), ( 1, miriam[ "synonym" ] ), ( 2, miriam[ "cas" ] ), ( 4, miriam[ "kegg.compound" ] ),
                                ( 5, miriam[ "kegg.drug" ] ), ( 6, miriam[ "pubchem.compound" ] ), ( 7, miriam[ "pubchem.substance" ] ),
                                ( 8, miriam[ "chebi" ] ), ( 9, miriam[ "pharmgkb.drug" ] ), ( 11, miriam[ "uniprot" ] )
                            )
                stream.readline( )  # skip header
                for record in stream:
                    line = record.split( sep )
                    try:
                        # column 12 contains the Uniprot Name (e.g. CDKN1_HUMAN)
                        try:
                            gene_id = uniprot[ line[ 11 ] ]
                        except KeyError:
                            gene_id = None
                        if not gene_id:
                            i += 1
                            gene_id = i
                            drugs[ gene_id ] = dict( symbol = line[ 1 ].strip( ), tax_id = 1, desc = "from Drugbank" )
                            xrefs[ gene_id ] = set( )
                        for p, ref_type in alias_pos:
                            if line[ p ].strip( ):
                                xrefs[ gene_id ].add( ( ref_type, line[ p ].strip( ) ), )
                    except IndexError:
                        self._alert( "Notice: Incompatible line in {!r}: {}" . format( db_file, line[ 0 ] ) )
                stream.close( )
            else:
                self._spill( "Unrecognized file name schema: {!r}" . format( db_file ) )

        # map miRNA aliases to the corresponding gene
        count = 0
        try:
            for mB_acc, al in aliased_mirnas.iteritems( ):
                try:
                    gene_id = mapped_mirnas[ mB_acc ]
                except KeyError:
                    count += 1
                    #self._spill( "Warning: no matching gene for miRBase accession {}." . format( mB_acc ) )
                    continue
                try:
                    if mB_acc[ :5 ] == "MIMAT":
                        xrefs[ gene_id ].add( ( miriam[ "mirbase.mature" ], mB_acc ), )
                    xrefs[ gene_id ] |= al
                except KeyError:
                    continue
            self._spill( "Info: {} miRBase identifiers lack gene annotation." . format( count ) )
        except NameError:  # no variable named mapped_mirnas because no miRBase data file was specified
            self._alert( "A miRBase alias file needs to be read together with the corresponding miRBase data file." )

        # update annotation table with newly added types
        if new_miriams:
            query = "INSERT INTO `Actor_xrefs` VALUES ({})" . format( "),\n\t(" . join( [ self._sqllist( m ) for m in new_miriams ] ) )
            self._sql( query )

        # prepare entries for role table
        roles = sorted( roles )
        role_dict = dict( zip( roles, list( xrange( 2, len( roles ) + 2 ) ) ) )
        insert_values[ tRoles ] = [ self.unknown_entity, "MISSING_VALUE" ]
        for e, v in role_dict.items( ):
            insert_values[ tRoles ].extend( [ v, e ] )

        # prepare gene entries for actors table
        i = 0
        entrez_to_id = { }
        for k, v in annot.iteritems( ):
            i += 1
            entrez_to_id[ k ] = i
            try:
                insert_values[ tAnnot ].extend( [ v[ "tax_id" ], i, v[ "symbol" ], v[ "desc" ] ] )
            except:
                print( tAnnot, k, v )
                sys.exit( )
            try:
                for ref in xrefs[ k ]:
                    insert_values[ tAlias ].extend( [ i, ref[ 1 ], ref[ 0 ] ] )
            except ( IndexError, KeyError ):
                self._spill( "Error: unmapped identifier: {}  {}  {}" . format( i, k, ref ) )
                sys.exit( )

        # the unknown entity needs to be treated separately
        insert_values[ tAnnot ].extend( [ 1, no_gene, "UNMAPPABLE_IDENTIFIER", "marks agents whose identity remains unclear" ] )
        entrez_to_id[ no_gene ] = no_gene
        for ref in xrefs[ no_gene ]:
            insert_values[ tAlias ].extend( [ no_gene, ref[ 1 ], ref[ 0 ] ] )

        # prepare complex entries for actor table
        for k, v in complexes.iteritems( ):
            insert_values[ tAnnot ].extend( [ v[ "tax_id" ], k, v[ "symbol" ], v[ "desc" ] ] )
            try:
                for ref in xrefs[ k ]:
                    insert_values[ tAlias ].extend( [ k, ref[ 1 ], ref[ 0 ] ] )
            except ( IndexError, KeyError ):
                self._spill( "Error: unmapped identifier: {}  {} " . format( k, ref ) )
                sys.exit( )
            for m, c in v[ "members" ].iteritems( ):
                insert_values[ tCompl ].extend( [ k, entrez_to_id[ m[ 0 ] ], role_dict[ m[ 1 ] ], c ] )

        # prepare drug entries for actor table
        for k, v in drugs.iteritems( ):
            try:
                insert_values[ tAnnot ].extend( [ v[ "tax_id" ], k, v[ "symbol" ], v[ "desc" ] ] )
            except:
                print( tAnnot, k, v )
                sys.exit( )
            try:
                for ref in xrefs[ k ]:
                    insert_values[ tAlias ].extend( [ k, ref[ 1 ], ref[ 0 ] ] )
            except ( IndexError, KeyError ):
                self._spill( "Error: unmapped identifier: {}  {}  {}" . format( i, k, ref ) )
                sys.exit( )

        # prepare miRBase table entries
        insert_values[ "mirbase" ] = [ ]
        for mi, miv in mirna_refs.iteritems( ):
            try:
                gene_id = entrez_to_id[ mapped_mirnas[ mi ] ]
            except KeyError:
                gene_id = entrez_to_id[ no_gene ]
            for mat in miv[ "derived" ]:
                insert_values[ "mirbase" ].extend( [ gene_id, gene_id, "," . join( miv[ "PMIDs" ] ), mi, mat ] )

        return insert_values


    ############################################################
    #### Drugbank                                           ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_Drugbank_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            stream.readline( )  # skip header
            for record in stream:
                line = record.split( sep )
                try:
                    drugs = line[ 12 ].strip( ).split( "; " )  # NOTE: needs translation
                    target = line[ 0 ].strip( )  # NOTE: needs translation

                    for drug in drugs:
                        interactions.add( ( drug, target, "" ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "drugbank" ], [ "uniprot" ] ) )

        return {  self.getTableConf( "Drugbank" )[ "name" ] : insert_values  }


    ############################################################
    #### HPRD                                               ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_HPRD_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            for record in stream:
                line = record.split( sep )
                try:
                    refseq1 = line[ 2 ].strip( ).split( "." )[ 0 ] # NOTE: needs translation
                    refseq2 = line[ 5 ].strip( ).split( "." )[ 0 ]  # NOTE: needs translation
                    validation = line[ 6 ].strip( )
                    pmids = line[ 7 ].strip( )

                    interactions.add( ( refseq1, refseq2, pmids ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "refseq" ], [ "refseq" ] ) )

        return {  self.getTableConf( "HPRD" )[ "name" ] : insert_values  }


    ############################################################
    #### HTRIdb                                             ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_HTRIdb_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror + ": " + db_file )
                return 1

            next( stream )  # discard headers

            for record in stream:
                line = record.split( sep )
                try:
                    id1 = line[ 1 ].strip( )  # NOTE: needs translation
                    id2 = line[ 3 ].strip( )  # NOTE: needs translation
                    pmids = line[ 6 ].strip( )

                    interactions.add( ( id1, id2, pmids ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "ncbigene" ], [ "ncbigene" ] ) )

        return {  self.getTableConf( "HTRIdb" )[ "name" ] : insert_values  }


    ############################################################
    #### miRecords                                          ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_miRecords_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            next( stream )  # discard file headers
            for record in stream:
                line = record.split( sep )
                try:
                    # in version 4, each line is comprised of 26 columns
                    pmid = int( line[ 0 ] )
                    #prot_name = line[ 2 ].strip( )  # NOTE: needs translation
                    refseq = line[ 3 ].strip( ).split( "." )[ 0 ]  # NOTE: needs translation
                    try:
                        site_num = int( line[ 4 ] )
                    except ValueError:
                        site_num = None
                    miRNA_name = line[ 6 ].strip( " \t[]" )  # NOTE: needs translation (caution: version 4 uses the deprecated asterisk nomenclature)
                    if miRNA_name[ :3 ] in ( "mir", "miR", "let" ):
                        species = get_latin_abbr( line[ 5 ] )
                        miRNA_name = species + "-" + miRNA_name
                    elif miRNA_name == "has-let-7a3b":
                        miRNA_name = "hsa-let-7a-5p"
                    elif miRNA_name == "has-let-71f1":
                        miRNA_name = "hsa-let-7f-5p"
                    elif miRNA_name[ :4 ] == "has-":
                        miRNA_name = "hsa-" + miRNA_name[ 4: ]
                    reg = line[ 7 ].strip( )  # regulation
                    rep_link_el = line[ 9 ].strip( )  # reporter_link_element
                    try:
                        site_pos = int( line[ 16 ] )
                    except ValueError:
                        site_pos = None
                    if reg in ( "NA", "" ):  reg = "intranslation"
                    if rep_link_el in ( "NA", "" ):  rep_link_el = None

                    interactions.add( ( miRNA_name, refseq, pmid ) )
                except ( ValueError, IndexError ):
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "synonym" ], [ "refseq" ] ) )

        return {  self.getTableConf( "miRecords" )[ "name" ] : insert_values  }


    ############################################################
    #### miRTarBase                                         ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_miRTarBase_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror + ": " + db_file )
                return 1

            next( stream )  # discard file header
            for record in stream:
                line = record.split( sep )
                try:
                    #acc = line[ 0 ].strip( )
                    mirna = line[ 1 ].strip( )  # NOTE: needs translation
                    gene_entrez = line[ 4 ].strip( )  # NOTE: needs translation
                    if not gene_entrez:
                        gene_entrez = line[ 3 ].strip( )  # use symbol instead
                    #proof = line[ 6 ].strip( )
                    #supp = line[ 7 ].strip( )
                    pmids = line[ 8 ].strip( )

                    # the 'species' column refers to the messenger-RNA's organism (e.g. the host in case of viruses), so 'gene_spec' WAS used
                    interactions.add( ( mirna, gene_entrez, pmids ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "synonym" ], [ "ncbigene", "hgnc.symbol" ] ) )

        #return {  self.getTableConf( "miRTarBase" )[ "name" ] : insert_values  }
        return {  self.getTableConf( "miRTarBase6" )[ "name" ] : insert_values  }


    ############################################################
    #### RegPhos                                            ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_RegPhos_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            next( stream )  # discard file header
            for record in stream:
                line = record.split( sep )
                try:
                    kinase = line[ 4 ].strip( )
                    pmids = "," . join( line[ 5 ].strip( ).split( ";" ) )
                    if not kinase or not pmids:
                        continue  # skip entries where no kinase or no reference is specified
                    symbol, species = line[ 0 ].strip( ).split( "_" )
                    uniprot_acc = line[ 1 ].strip( )
                    #phosphoposition = int( line[ 2 ].strip( ) )
                    #residue = line[ 3 ].strip( )[ 7:10 ]

                    if symbol.strip( ) == "autocatalysis":
                        symbol = uniprot_acc

                    interactions.add( ( kinase, uniprot_acc, pmids ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "hgnc.symbol", "uniprot" ], [ "uniprot" ] ) )

        return {  self.getTableConf( "RegPhos" )[ "name" ] : insert_values  }


    ############################################################
    #### starBase                                           ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_starBase_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        insert_values = [ ]
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            next( stream )  # discard file header
            for record in stream:
                line = record.split( sep )
                try:
                    mirna = line[ 0 ].strip( )  # NOTE: needs translation
                    gene = line[ 1 ].strip( )  # NOTE: needs translation
                    targetScan = line[ 2 ].strip( )
                    picTar = line[ 3 ].strip( )
                    RNA22 = line[ 4 ].strip( )
                    PITA = line[ 5 ].strip( )
                    miRanda = line[ 6 ].strip( )

                    insert_values.extend( [ mirna, gene, targetScan, picTar, RNA22, PITA, miRanda ] )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )
        return {  self.getTableConf( "starBase" )[ "name" ] : insert_values  }


    ############################################################
    #### STRING                                             ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_STRING_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        insert_values = [ ]
        return [ ]  # TODO: remove after implementation is done

        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            for record in stream:
                line = record.split( sep )
                try:

                    insert_values.extend( [  ] )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )
        return {  self.getTableConf( "STRING" )[ "name" ] : insert_values  }


    ############################################################
    #### TarBase                                            ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_TarBase_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            for record in stream:
                line = record.split( sep )
                try:
                    target = line[ 3 ].strip( ).strip( '"' )  # NOTE: needs translation
                    target = target.split( "/" )[ 0 ].strip( )  # remove second name in double names
                    #action_type = line[ 0 ].strip( )
                    pmid = line[ 4 ].strip( )
                    try:
                        int( pmid )
                    except ValueError:
                        pmid = 0
                    species = get_scientific_name( line[ 1 ] )  # NOTE: needs translation
                    mirna = "{}-{}" . format( get_latin_abbr( species ), line[ 2 ].strip( ) )  # NOTE: needs translation

                    interactions.add( ( mirna, target, pmid ) )
                except IndexError:
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "synonym" ], [ "hgnc.symbol" ] ) )

        return {  self.getTableConf( "TarBase" )[ "name" ] : insert_values  }


    ############################################################
    #### TRANSFAC                                           ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_TRANSFAC_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            next( stream )  # discard file header
            for record in stream:
                line = record.split( sep )
                try:
                    species = line[ 3 ].strip( )  # NOTE: needs translation
                    if species == "Homo sapiens":
                        species = 9606
                    elif species == "Mus musculus":
                        species = 10090
                    else:
                        continue
                    gene_id = line[ 2 ].strip( )  # NOTE: needs translation
                    factor_ids = line[ 5 ].strip( )  # NOTE: needs translation
                    pmids = line[ 6 ].strip( )

                    for f in factor_ids.split( "," ):
                        interactions.add( ( f, gene_id, pmids ) )
                except ( IndexError, ValueError ):
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "transfac" ], [ "transfac" ] ) )

        return {  self.getTableConf( "TRANSFAC" )[ "name" ] : insert_values  }


    ############################################################
    #### TransmiR                                           ####
    ############################################################

    ### encapsulated method for reading data from a text file ###
    def _read_TransmiR_file( self, db_files, sep ):
        """read the database file and return a dict of VALUE lists for the INSERT query"""

        interactions = set( )
        for db_file in db_files:
            try:
                stream = open( db_file, "r" )
            except IOError as e:
                self._extalert( e.strerror )
                return 1

            next( stream )  # discard file header
            for record in stream:
                line = record.split( sep )
                try:
                    gene_id = line[ 1 ].strip( )  # NOTE: needs translation
                    #action_type = line[ 7 ].strip( ).lower( )
                    pmids = "," . join( line[ 8 ].strip( ).split( ";" ) )
                    mirna = "{}-{}" . format( get_latin_abbr( line[ 9 ] ), line[ 3 ].strip( ) )  # NOTE: needs translation

                    interactions.add( ( gene_id, mirna, pmids ) )
                except ( IndexError, ValueError ):
                    self._alert( "Notice: Incompatible line in {!r}." . format( db_file ) )

            stream.close( )

        insert_values = self._resolve_identifiers( interactions, ( [ "ncbigene" ], [ "synonym" ] ) )

        return {  self.getTableConf( "TransmiR" )[ "name" ] : insert_values  }


############################################################
#### QUERY METHODS                                      ####
############################################################


    def _resolve_identifiers( self, interactions, types ):
        """query the database and translate identifiers before writing them to interaction tables"""

        self._spill( "\tResolving identifiers..." )
        srcs, tgts = zip( *interactions )[ :2 ]
        s = self.getAllowedSpecies( )
        self.setAllowedSpecies( None )
        srcs_tr = self.query_for_references( srcs, { "alias_types" : types[ 0 ] }, invert = True )
        tgts_tr = self.query_for_references( tgts, { "alias_types" : types[ 1 ] }, invert = True )
        self.setAllowedSpecies( s )

        insert_values = [ ]
        for s, t, p in interactions:
            # prepare shared identifiers
            try:  # success one-to-one?
                si = srcs_tr[ 0 ][ s.lower( ) ],
            except KeyError:
                try:  # success one-to-many?
                    si = srcs_tr[ 2 ][ s.lower( ) ]
                except KeyError:  # orphan
                    si = self.unknown_entity,
            try:
                ti = tgts_tr[ 0 ][ t.lower( ) ],
            except KeyError:
                try:
                    ti = tgts_tr[ 2 ][ t.lower( ) ]
                except KeyError:
                    ti = self.unknown_entity,
            for source in si:
                for target in ti:
                    insert_values.extend( [ source, target, p, s, t ] )

        return insert_values


    ############################################################
    #### PUBLIC                                             ####
    ############################################################

    def setAllowedSpecies( self, species_list ):
        """configure and return the set of species used for filtering subsequent SQL queries"""

        if not species_list:
            self._allowed_species = [ ]
            return
        elif type( species_list ) in ( list, tuple, dict, set ):
            self._allowed_species = list( species_list )
        elif type( species_list ) in ( str, int ):
            self._allowed_species = [ species_list ]
        else:
            self._extalert( "Parameter species_list: expected iterable, set, str, int, or None." )
            return 1

        select = [ "`ref`" ]
        table = {  "`Taxon_aliases`" : ""  }
        cond = [  {  "`Alias`" : [ "IN" ] + [  str( s ) for s in self._allowed_species  ] }  ]
        try:
            self._allowed_species = set( *self._select( select, table, cond ) )
        except TypeError:
            self._alert( "None of the specified species aliases was found in the database." )
            return 2

        return self._allowed_species


    def getAllowedSpecies( self ):
        """return the set of species used for filtering subsequent SQL queries"""

        return self._allowed_species


    ############################################################
    #### The query skeletons below are structured like this:
    #### 1. Annotation queries
    ####    - single key + dictionary 2-tuples = uneven column count
    #### 2. Interaction queries
    ####    - two keys + two types + dictionary 2-tuples = even column count
    ############################################################

    def _extract_dict( self, row, fixed ):
        """interpret all elements from position 'fixed' on as two-tuples for a dictionary"""

        r = row
        s = fixed
        return dict( [ ( r[ s + 2 * n ], r[ s + 2 * n + 1 ] ) for n in xrange( ( len( r ) - s ) / 2 ) ] )


    ############################################################
    #### aliases                                            ####
    ############################################################

    def query_for_references( self, candidates, restrict = { "alias_types" : None }, invert = False, silent = True ):
        """retrieve shared identifiers for the given aliases"""

        identified = { }  # will hold candidates that mapped to exactly one shared identifier (one-to-one)
        ambiguous = { }  # will hold candidates that mapped to multiple shared identifers (one-to-many)
        orphans = set( )  # will hold candidates that were unmappable (one-to-zero)
        if type( candidates ) not in ( tuple, list, dict, set ):
            self._exalert( "Wrong parameter syntax, aborting: expected iterable for parameter 'candidates'" )
            return 1
        if type( restrict ) != dict or "alias_types" not in restrict:
            restrict = { "alias_types" : None }
        elif type( restrict[ "alias_types" ] ) == str:
            restrict[ "alias_type" ] = set( [ restrict[ "alias_types" ] ] )

        requested = {  str( e ).lower( ) for e in candidates  }
        if not requested:
            return ( identified, orphans, ambiguous )  # return empty dict for empty candidate list

        c = 0
        db = "Actor_aliases"
        if invert:
            order = "`Alias`, `ref`"
            p = { "ref" : 1, "Alias" : 0 }
        else:
            order = "`ref`, `Alias`"
            p = { "ref" : 0, "Alias" : 1 }
        query = "SELECT\n\t{}\nFROM\n\t`Actor_aliases`\nWHERE\n\t" . format( order )
        if restrict[ "alias_types" ]:
            query += "`type` IN ( SELECT `x_id` FROM `Actor_xrefs` WHERE `namespace` IN ( {} ) )\n\tAND\n\t" . format( ", " . join( self._sqlescape( restrict[ "alias_types" ] ) ) )
        query += "`Alias` IN ( {} )" . format( ", " . join( self._sqlescape( list( requested ) ) ) )
        if self._allowed_species:
            query = "SELECT\n\t{}\nFROM\n\t( {} ) AS `a`\n\tJOIN `Actors` AS `g`\n\tON `a`.`ref` = `g`.`a_id`\nWHERE\n\t`g`.`species` IN ( {} )" . format( order, query.replace( "\n", "\n\t" ), ", " . join( self._sqlescape( self._allowed_species ) ) )
        res = self._select( query )

        if res == None:
            self._alert( "Misformatted SQL query, exiting." )
            sys.exit( 1 )

        if len( res ) == 0:
            orphans = requested  # nothing was found -> the whole input is classified as orphaned
        else:
            for alias in candidates:
                lalias = alias.lower( )
                s = {  t[ p[ "ref" ] ] for t in res if t[ p[ "Alias" ] ].lower( ) == lalias  }
                if len( s ) > 1:  # spreading - one alias maps to more than one reference (can happen for histone 3, for example)
                    try:
                        ambiguous[ alias ] |= s
                    except KeyError:
                        ambiguous[ alias ] =  s
                elif len( s ) == 0:  # orphans - aliases without a matching reference
                    orphans.add( alias )

            if not silent and ambiguous:
                self._spill( "Aliases invoking more than one reference: {:d} (out of {:d})." . format( len( ambiguous ), len( candidates ) ) )
            if not silent and orphans:
                self._spill( "Unresolvable aliases under given restrictions: {:d} (out of {:d})." . format( len( orphans ), len( candidates ) ) )

            # funneling - two or more aliases map to the same reference
            l = len( set( zip( *res )[ p[ "ref" ] ] ) )
            if l < len( res ):
                if not silent:
                    self._spill( "Identifier funneling reduced expected references by {:d}." .format( len( res ) - l ) )
                positions = { }
                for r in res:
                    try:
                        if r[ p[ "Alias" ] ].lower( ) not in {  a.lower( ) for a in positions[ r[ p[ "ref" ] ] ]  }:
                            positions[ r[ p[ "ref" ] ] ].append( r[ p[ "Alias" ] ] )
                    except KeyError:
                        positions[ r[ p[ "ref" ] ] ] = [ r[ p[ "Alias" ] ] ]
                if invert:
                    for k, v in positions.items( ):
                        identified.update( [ ( a.lower( ), k ) for a in v ] )
                else:
                    for k, v in positions.items( ):
                        identified[ k ] = " " . join( v )
            else:
                if invert:
                    identified = dict( [  ( e1.lower( ), e2 ) for e1, e2 in res  ] )
                else:
                    identified = dict( res )

        return ( identified, orphans, ambiguous )


    def query_for_aliases( self, references, restrict = None, silent = False ):
        """retrieve all aliases for the given shared identifiers"""

        aliased = { }
        if type( references ) not in ( tuple, list, dict, set ):
            self._extalert( "Wrong parameter syntax, aborting: expected iterable for parameter 'references'" )
            return 1

        requested = set( references )
        if not requested:
            return aliased

        query = """SELECT
    ref, type, Alias
FROM
    Actor_aliases
WHERE
    ref IN ( {} )"""

        if restrict:
            if type( restrict ) not in ( tuple, list, dict, set ):
                self._extalert( "Wrong parameter syntax, aborting: expected iterable for parameter 'restrict'" )
                return 1

            xrefs = dict( self._select( "SELECT `x_id`, `namespace` FROM `Actor_xrefs` WHERE `namespace` IN ( {} )" . format( ", " . join( self._sqlescape( restrict ) ) ) ) )
            if len( xrefs ) == len( set( restrict ) ):
                query += " AND type IN ( {} )" . format( ", " . join( self._sqlescape( xrefs.keys( ) ) ) )
            else:
                self._alert( "Warning: some xref types were funneled or not found: " + " " . join( restrict ) )

        res = self._select( query . format( ", " . join( self._sqlescape( requested ) ) ), None, None )

        d = aliased  # just another name to make typing easier
        found = set( )
        for r in res:
            if not r[ 1 ] in xrefs:
                continue
            try:
                d[ r[ 0 ] ][ xrefs[ r[ 1 ] ] ].append( r[ 2 ] )
            except KeyError:
                try:
                    d[ r[ 0 ] ][ xrefs[ r[ 1 ] ] ] = [ r[ 2 ] ]
                except KeyError:
                    d[ r[ 0 ] ] = {  xrefs[ r[ 1 ] ] : [ r[ 2 ] ]  }
            found.add( r[ 0 ] )

        # keep track of unmatched references
        if len( requested ) > len( found ):
            orphans = requested - found
            if not silent:
                self._spill( "No aliases with these restrictions for {:d} reference(s) (out of {:d})." . format( len( orphans ), len( requested ) ) )
        else:
            orphans = set( )


        return ( aliased, orphans )


    ############################################################
    #### annotations                                        ####
    ############################################################

    def query_for_annotations( self, candidates, silent = False ):
        """retrieve annotation data for the given identifiers"""

        helper = { }
        if type( candidates ) not in ( tuple, list, dict, set ):
            self._exalert( "Wrong parameter syntax, returning empty result." )
            return helper

        id_list = candidates
        db = "Actors"
        t = self.getTableConf( db )
        try:
            prefix = t[ "name" ]
        except ( KeyError, TypeError ):
            self._extalert( "Warning: Table configuration for {!r} not found." . format( db ) )
            #continue
            return helper

        selects = [
                "`a_id`",  # entry id (fixed)
                ]
        fixed = len( selects )
        selects += [
                "'symbol'", "`symbol`",
                "'description'", "`description`",
                "'species'", "`species`",
                ]
        #selects.append( "`{" + prefix + "_i002}_db`" )
        #for i in t[ "retrieve" ]:
        #    i = "{:03d}" . format( i )
        #    selects.extend( [ "'{" + prefix + "_i" + i + "}'", "`{" + prefix + "_i" + i + "}`" ] )  # name and value for additional columns

        froms = {  db : ""  }
        conds = [  {  "`a_id`" : [ "IN" ] + list( id_list )  }  ]

        res = self._select( selects, froms, conds )

        lowered_list = id_list
        found_list = list( res )
        queried = candidates
        found = set( found_list )
        # TODO:
        if not silent and len( found ) < len( queried ):
            self._spill( "Unable to annotate {:d} identifier(s) (out of {:d})." . format( len( queried ) - len( found ), len( queried ) ) )

        #helper[ mol_type ] = { }
        for r in res:
            helper[ r[ 0 ] ] =  self._extract_dict( r, 1 )

        return helper


    ############################################################
    #### interactions                                       ####
    ############################################################

    def query_for_interactions( self, candidates, restrict = {  "hooks" : None, "fish" : None, "ints" : None, "dbs" : None, "go" : None, "role" : None  } ):
        """retrieve interactions for the given identifiers, according to set of restrictions"""

        # GO terms database restriction:
        # 1. forwards:
        #SELECT `name`, `term_type`, `acc` FROM `term` WHERE id IN ( SELECT DISTINCT `term_id` FROM `association` WHERE `is_not` = 0 AND `gene_product_id` IN ( SELECT `id` FROM `gene_product` WHERE `species_id` = 355777 AND `id` IN ( SELECT `gene_product_id` FROM `gene_product_synonym` WHERE `product_synonym` IN ( 'VEGF' ) ) ) )
        # 2. backwards:
        #SELECT xref_dbname, xref_key FROM `dbxref` WHERE id IN ( SELECT dbxref_id FROM `gene_product` WHERE species_id = 355777 AND id IN ( SELECT DISTINCT gene_product_id FROM `association` WHERE is_not = 0 AND term_id IN ( SELECT id FROM `term` WHERE acc IN ( 'GO:0006955' ) ) ) )

        helper = [ ]
        if type( candidates ) not in ( tuple, list, dict, set ):
            self._exalert( "Wrong parameter syntax, returning empty result." )
            return helper
        elif not candidates:
            return ( )

        dbi_defs = self.getConfItem( "db_types" )[ "interactions" ]  # maps mol_type to src-tgt relation with databases
        int_defs = self.getConfItem( "interdefs" )  # maps interactions to src-tgt definitions
        query_defs = self._db_restriction  # allowed databases
        role_defs = "src", "tgt",  "two-way"
        for key, default in ( "hooks", dbi_defs ), ( "fish", dbi_defs ), ( "ints", int_defs ), ( "dbs", query_defs ), ( "go", [ ] ), ( "role", role_defs ):
            if key not in restrict or restrict[ key ] == None:
                restrict[ key ] = set( default )
            elif type( restrict[ key ] ) != str:
                try:
                    restrict[ key ] = set( restrict[ key ] )  # works only if key element is an iterable
                except:
                    self._alert( "The {!r} restriction is not a string and not iterable either, aborting." . format( key ) )
                    return helper
            else:
                restrict[ key ] = set( [ restrict[ key ] ] )  # string becomes the set's single element

        _org = self.getConfItem( "suffixes" )[ "org" ]
        id_list = list( set( candidates ) - set( [ self.unknown_entity ] ) )
        for mol_type in restrict[ "hooks" ]:  # only nodes of eligible types may serve as hooks
            for role, oppos, pos in ( "src", "tgt", ( "001", ) ), ( "tgt", "src", ( "002", ) ), ( "two-way", None, ( "001", "002" ) ):
                if role not in dbi_defs[ mol_type ] or role not in restrict[ "role" ]:
                    continue
                for db in dbi_defs[ mol_type ][ role ]:
                    t = self.getTableConf( db )
                    if (  # filter databases by...
                        t[ "name" ] not in restrict[ "dbs" ]  # ... name
                        or oppos != None and t[ oppos ] not in restrict[ "fish" ]   # ... mol_type to be added
                        or t[ "src" ] not in set.union( *[ int_defs[ s ][ "src" ] for s in restrict[ "ints" ] ] )  # ... interaction type (partner 1)
                        or t[ "tgt" ] not in set.union( *[ int_defs[ s ][ "tgt" ] for s in restrict[ "ints" ] ] )  # ... interaction type (partner 2)
                    ):
                        continue  # skip this table because it is not eligible for querying
                    try:
                        prefix = t[ "affix" ]
                    except KeyError:
                        self._extalert( "Warning: Table configuration for {!r} not found." . format( db ) )
                        continue
                    selects = [
                            "`source`", # source id (fixed)
                            "`target`", # target id (fixed)
                            "'" + t[ "src" ] + "'",  # source type (fixed)
                            "'" + t[ "tgt" ] + "'",  # target type (fixed)
                            ]
                    fixed = len( selects )  # the number of leading non-annotation columns in the query result (comes from definition of 'selects' )
                    selects = [ "DISTINCT" ] + selects + [
                            "'source_alias'", "`source_orig`",
                            "'target_alias'", "`target_orig`",
                            "'database'", "'" + db + "'",
                            "'release'", "'" + t[ "release" ] + "'",
                            "'reldate'", "'" + t[ "reldate" ] + "'",
                            "'PMIDs'", "`PMIDs`",
                            ]
                    #for i in t[ "retrieve" ]:
                    #    i = "{:03d}" . format( i )
                    #    selects.extend( [ "'{" + prefix + "_i" + i + "}'", "`{" + prefix + "_i" + i + "}`" ] )  # name and value for additional columns
                    froms = {  db : ""  }
                    conds = [ {  "{" + prefix + "_i" + p + "}" : [ "IN" ] + id_list  } for p in pos ]
                    res = self._select( selects, froms, conds )
                    for r in res:
                        helper.append( list( r[ :fixed ] ) + [ self._extract_dict( r, fixed ) ] )
                    #query = "SELECT\n\tsource, target, {}, {}, 'source_alias', source_orig, 'target_alias', target_orig, 'database"

        return helper


    ############################################################
    #### GENERAL AND NOT SO PUBLIC                          ####
    ############################################################

    def _sqlescape( self, elements ):
        """escape the elements for an sql query"""
        return self._mysql_escape( list( elements ), self._mysql_escape_conversions )

    def _sqllist( self, elements ):
        """sql-escape elements and concatenate into string"""
        return ", " . join( self._mysql_escape( list( elements ), self._mysql_escape_conversions ) )


    def _select( self, sql_sel, sql_from = None, sql_where = None ):
        """handle arbitrary SELECT statements (by parsing a complicated syntax)"""

        if [ type( sql_sel ), type( sql_from ), type( sql_where ) ] == [ list, dict, list ]:

            # The data structures are parsed according to the rules outlined below.
            # Alias specifications are always optional and mainly good for checking complex queries by eye.

            # build the SELECT clause
            # sql_sel list elements are strings '{column_identifier}:alias', i.e. "SELECT column_name1 AS alias1[, column_name2 AS alias2[, ...]]"
            if sql_sel[ 0 ] == "DISTINCT":
                sel_part = "SELECT DISTINCT\n\t"
                del sql_sel[ 0 ]
            else:
                sel_part = "SELECT\n\t"
            helper = [ ]
            for s in sql_sel:
                fields = s.split( ":" )
                if len( fields ) > 1:
                    helper.append( "{} AS {}" . format( fields[ 0 ], fields[ 1 ] ) )
                else:
                    helper.append( "{}" . format( fields[ 0 ] ) )
            sel_part += ",\n\t" . join( helper )

            # build the FROM clause
            # sql_from dictionary contains dictionaries { '{table_identifier}' : 'alias' }, i.e. "FROM table1[ AS alias1][, table2[ AS alias2[, ...]]]"
            helper = [ ]
            for k, v in sql_from.items( ):
                if k == "JOIN":
                    ( t1, t1a ) = v[ "T1" ].items( ).pop( )
                    ( t2, t2a ) = v[ "T2" ].items( ).pop( )
                    if t1a:
                        t1a = " AS " + t1a
                    if t2a:
                        t2a = " AS " + t2a
                    helper.append( "{}{} JOIN {}{} ON {}" . format( t1, t1a, t2, t2a, v[ "ON" ] ) )
                elif v != "":
                    helper.append( "{} AS {}" . format( k, v ) )
                else:
                    helper.append( "{}" . format( k ) )
            from_part = "FROM\n\t" + ",\n\t" . join( helper )

            # build the WHERE clause
            # sql_where list elements are the most complex
            # - they are themselves dictionaries with list values: list of dict of list == [ { 'key1' : [ value1, value2, ... ], 'key2' : ... }, { ... } ]
            # - each dictionary is interpreted as { 'column_ref' : [ 'comparison_mode', value1[, value2[, ...]] }, all comparisons are disjunct (OR)
            # - 'comparison_mode' is one of '=', '!s', 'LIKE', 'NOT LIKE' (for these, the other elements will be quoted) or any other valid MySQL comparator
            # - the outer list elements are then connected to form a single conjunction (AND) (now it forms a union - need to set up a viable logic for this)
            conjunctions = [ ]
            for cond in sql_where:
                helper = [ ]
                for k, v in cond.items( ):
                    if len( v ) == 1:
                        continue
                    try:
                        v[ 1: ] = self._mysql_escape( v[ 1: ], self._mysql_escape_conversions )
                        if v[ 0 ] in ( "IN", "NOT IN" ):
                            helper.append( "{} {} ( {} )" . format( k, v[ 0 ], ", " . join( v[ 1: ] ) ) )
                        else:
                            helper.extend( [ "{} {} {}" . format( k, v[ 0 ], subj ) for subj in v[ 1: ] ] )
                    except ( TypeError, IndexError ):
                        self._extalert( "Skipping condition entry {!r} in query dictionary." . format( k ) )
                conjunctions.append( "( " + " OR " . join( helper ) + " )" )

            where_part = "WHERE\n\t" + "\n\tOR " . join( conjunctions )  # NOTE: changed AND to OR

            # combine all three and format the query
            query = str( "\n" . join( [ sel_part, from_part, where_part ] ) ) . format( **self.getSQLMapping( ) )

        elif type( sql_sel ) == str and sql_from == None and sql_where == None:
            query = sql_sel . format( **self.getSQLMapping( ) )
        else:
            self._alert( "Warning: Method '_select' takes either list, dict, list or one character string as arguments." )
            return
        #print( query, "\n" )

        try:
            self._cursor.execute( query )
            rows = self._cursor.fetchall( )
            #print( "Queried {}, {:d} hits" . format( ", ".  join( sql_from.keys( ) ), len( rows ) ) )
            return rows
        except MySQLdb.Error as e:
            self._alert( "SQL Query Error.\nDescription: {}\nAffected SQL Query:\n{}\n\n" . format( e, query ) )
            return ( )


    def _sql( self, sql_statement ):
        """execute arbitrary sql statements"""

        try:
            self._cursor.execute( sql_statement )
            return self._cursor.fetchall( )  # collect all rows that the query yielded
        except MySQLdb.Error as e:
            self._alert( "SQL Error.\nDescription: {}\nAffected SQL Query:\n{}\n\n" . format( e, sql_statement ) )
            return 1

