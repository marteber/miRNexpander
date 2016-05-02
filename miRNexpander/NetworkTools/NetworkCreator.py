#!/usr/bin/env python

from __future__ import print_function

import MySQLdb
import os, sys
import re
import networkx
from collections import defaultdict  # for easier handling of dictionaries with nested entries
from copy import deepcopy  # for copying nested types by value
from itertools import chain  # for selecting any element from a group of iterables
#from networkx import gviz_layout
#import matplotlib.pyplot as plt
from networkx.readwrite import json_graph  # for json export
import json  # for json export
from miRNexpander.NetworkTools.AliasResolver import AliasResolver
from miRNexpander.DatabaseTools.DatabaseHandler import DatabaseHandler
from miRNexpander.SBMLTools.SBMLTools import NetworkWriter
from miRNexpander.QueryTemplates.QT import query_templ



class NetworkCreator( AliasResolver ):
    """build a network from given seeds"""

    # the following dictionary is a template; it maps HGNC symbols to nodes and expression data, using different categories
    __symbol_mapper = dict(
        node = defaultdict(  # gene level
            lambda: defaultdict( list ),  # node level (base, active, phospho, ...)
        ),
        expr = defaultdict(  # gene level
            lambda: defaultdict( float ),  # treatment group level (Lpne_cont, Lpne_inf, Mtub_cont, Mtub_inf, ...)
        )
    )
    __mirna_regexp = re.compile( "mir(let)?[0-9]+[a-z]?(-[0-9]+)?" )
    __hostgene_regexp = re.compile( "mir(let)?[0-9]+[a-z]?(-[0-9]+)?hg" )

    def __init__( self, db_handler ):
        """initialize instance variables"""

        self.graph = networkx.MultiDiGraph( )
        self.Name = None

        # declare the use copy of the class variable dictionary
        self._symbol_mapper = deepcopy( NetworkCreator.__symbol_mapper )

        self._db = db_handler
        self._moltypes = self._db.getConfItem( "moltypes" )
        self._output_handler = None
        self._reg_type = self._db._reg_type
        self._symmetric =  ( "binding", "association", "co-avtivation", "cross-activation", )


############################################################
#### PUBLIC                                             ####
############################################################

    def setName( self, name ):
        """set the network's name"""
        if type( name ) == str:
            self.Name = name
            self.graph.graph[ 'name' ] = name


    def createNetwork( self, seeds, depth = 1, restrict = { "fish" : None }, cfilter = "any" ):
        """query the database to build a network"""

        if not seeds or type( seeds ) not in ( tuple, list, dict, set ):
            self._extalert( "Parameter 1 is empty or not an iterable." )
            return False
        try:
            int( depth )
        except ValueError:
            self._extalert( "Parameter 2 is not an integer value." )
            return False

        self._reset( )
        self.graph.graph[ "name" ] = self.Name
        # ToDo: check logic of species restriction
        self._add( seeds, seeds = True, complex_filter = cfilter )
        self._connect( add_shells = depth, cfilter = cfilter, restrict = restrict )
        return True


    def create_filterlist( self, means = { "GO.leafs" } ):
        """create a gene list for filtering"""

        gfilter = dict( )

        # prepare species
        sp = self._db.getAllowedSpecies( )
        if sp:
            sp_switch = "0"  # False -> species subclause must be evaluated
            sp = self._db._sqllist( sp )
        else:
            sp_switch = "1"  # True -> species subclause can be skipped
            sp = "0"  # irrelevant because
        fillings = dict(
                GO = "GO_20150711",  # GO database name
                MX = "miRNexpander",  # miRNexpander database name
                sp_switch = sp_switch,
                species = sp,
                namespace = "'hgnc.symbol'",
                input = self._db._sqllist( self._symbol_mapper[ "node" ].keys( ) )
        )
        for fi in means:
            try:
                query = query_templ[ fi ]
            except KeyError:
                self._spill( "No such query template: {!r}." . format( fi ) )
                continue
            # customize query with current input
            query = query.format( **fillings )
            # MySQLdb must handle one query after the other, so run the preparation statements and ignore their output.
            queries = query.split( ";" )
            for q in queries[ :-1 ]:
                self._db._sql( q )
            # the final statement is supposed to yield the result
            gfilter[ fi ] = set( self._db._sql( queries[ -1 ] ) )

        return gfilter


    def writeCSV( self, path = None ):
        """generate and write a basic CSV format"""

        if not self._output_handler:
            self._output_handler = NetworkWriter( self.graph )
            if path == None:
                path = os.path.join( self._db.getConfItem( "output_path" ) )
            self._output_handler.setOutputDir( path )
        self._output_handler.createCytoscape( )


    def writeSBML( self, path = None ):
        """generate and write the SBML format"""

        if not self._output_handler:
            self._output_handler = NetworkWriter( self.graph )
            if path == None:
                path = os.path.join( self._db.getConfItem( "output_path" ) )
            self._output_handler.setOutputDir( path )
        self._output_handler.createSBML( )


    def writeGML( self, path = None ):
        """generate and write the GML format"""

        if not self._output_handler:
            self._output_handler = NetworkWriter( self.graph )
            if path == None:
                path = os.path.join( self._db.getConfItem( "output_path" ) )
            self._output_handler.setOutputDir( path )
        self._output_handler.createGML( gml_dict = self._db.getGraphicsConf( "xgmml" ) )


    def writeXGMML( self, path = None ):
        """generate and write the XGMML format"""

        if not self._output_handler:
            self._output_handler = NetworkWriter( self.graph )
            if path == None:
                path = os.path.join( self._db.getConfItem( "output_path" ) )
            self._output_handler.setOutputDir( path )
        self._output_handler.createXGMML( xgmml_dict = self._db.getConfItem( "graphics" )[ "xgmml" ] )


    def snapshot( self, ft = "png" ):
        """take a picture of the network"""

        if not self._output_handler:
            self._output_handler = NetworkWriter( self.graph )
            if path == None:
                path = os.path.join( self._db.getConfItem( "output_path" ) )
            self._output_handler.setOutputDir( path )
        self._output_handler.saveImage( )


    def readNetworkFiles( self, mirna_file_list, protein_file_list ):
        """reads a network from files and marks all nodes as seeds"""

        return  # NOTE: not tested in this version

        # internal subroutine
        def _read_files( file_list, mol_type, target_dict ):
            """reads nodes of a specific type from files"""

            mols = list( )
            for f in file_list:
                try:
                    stream = open( f )
                except IOError:
                    self._extalert( "Error reading {!r}." . format( f ) )
                    continue

                for line in stream:
                    mol = line.strip( )
                    if mol:
                        mols.append( mol )

            if mol_type == "miRNA":
                self._add_miRNA_nodes( list( set( mols ) ), target_dict )
            elif mol_type == "protein":
                self._add_protein_nodes( list( set( mols ) ), target_dict )
            if not result_set:
                self._alert( "SQL Miss: No matching database entries for {!r}." . format( "," . join( seeds ) ) )
                return 1

        _read_files( mirna_file_list, "miRNA", self.Nodes )
        _read_files( protein_file_list, "protein", self.Nodes )


    ### FILE OUTPUT
    def write_graph( self, path, f = [ "gml" ], graph = None ):
        """write graph to file in specified format"""

        if type( f ) not in ( tuple, list, set, dict ):
            self._alert( "Not an iterable: {}" . format( repr( f ) ) )
            return

        if graph == None:
            graph = self.graph

        if "gml" in f:  # Graph Markup Language
            fn = path + ".gml"
            networkx.write_gml( graph, fn )
        if "d3.json" in f:  # JavaScript Object Notation, in Cytoscape.js-compatible format
            fn = path + ".json"
            jsong = graph.copy( )
            for n in jsong.nodes_iter( ):
                jsong.node[ n ][ "x" ] = jsong.node[ n ][ "graphics" ][ "x" ]
                jsong.node[ n ][ "y" ] = jsong.node[ n ][ "graphics" ][ "y" ]
            try:
                stream = open( fn, "w" )
            except IOError:
                self._spill( "Unable to open {} for writing." . format( fn ) )
            else:
                json.dump( json_graph.node_link_data( jsong ), stream )
                stream.close( )
        if "cy.json" in f:  # JavaScript Object Notation, in D3.js-compatible format
            fn = path + ".cyjs"
            jsong = networkx.MultiDiGraph( )
            jsong.add_edges_from( self.graph.edges_iter( data = True ) )
            for n in jsong.nodes_iter( ):
                jsong.node[ n ] = dict( position = dict( x = graph.node[ n ][ "graphics" ][ "x" ], y = graph.node[ n ][ "graphics" ][ "y" ] ),
                                        #data = dict( shared_name = graph.node[ n ][ "tag" ], name = graph.node[ n ][ "name" ], SUID = graph.node[ n ][ "tag" ] ),
                                        data = graph.node[ n ],
                                        selected = False
                                        )
            try:
                stream = open( fn, "w" )
            except IOError:
                self._spill( "Unable to open {} for writing." . format( fn ) )
            else:
                json.dump( json_graph.node_link_data( jsong ), stream )
                stream.close( )
        if "edges.csv" in f:  # "source interaction target" 3-column format for Cytoscape imports
            fn = path + "_edges.csv"
            try:
                attrib = f[ "edges.csv" ]
            except ( KeyError, TypeError ):
                attrib = "polarity"
            try:
                stream = open( fn, "w" )
            except IOError:
                self._spill( "Unable to open {} for writing." . format( fn ) )
            else:
                stream.write( "Source\tInteraction\tTarget\n" . format( attrib ) )
                stream.writelines( [ "{}\t{}\t{}\n" . format( e[ 0 ], e[ 2 ][ attrib ], e[ 1 ] ) for e in graph.edges_iter( data = True ) ] )
                stream.close( )
        if "adj.matrix" in f:  # adjacency matrix (node x node, value describes influence of row on column)
            fn = path + "_adj.matrix"
            try:
                stream = open( fn, "w" )
            except IOError:
                self._spill( "Unable to open {} for writing." . format( fn ) )
            else:
                nodes = sorted( graph.nodes( ) )
                n_ind = dict( [ ( nodes[ i ], i ) for i in xrange( len( nodes ) ) ] )
                stream.write( "\t" . join( nodes ) + "\n" )
                for n in nodes:
                    l = [ "0" ] * len( nodes )
                    for target, instances in graph.edge[ n ].iteritems( ):
                        # we need to check if all instances have the same polarity
                        try:
                            s = set( [ inst[ "polarity" ] for inst in instances.values( ) ] )
                        except KeyError:
                            print( graph.edge[ n ][ target ] )
                        else:
                            if len( s ) == 1:  # polarities agree or only one instance -> use polarity
                                l[ n_ind[ target ] ] = str( s.pop( ) )
                            else:  # polarities disagree
                                l[ n_ind[ target ] ] = "NA"
                    stream.write( "\t" . join( [ n ] + l ) + "\n" )
                stream.close( )
        if "fba.matrix" in f:
            fn = path + "_fba.matrix"
            try:
                stream = open( fn, "w" )
            except IOError:
                self._spill( "Unable to open {} for writing." . format( fn ) )
            else:
                # prepare dictionary to hold accessory representatives for homomultimeric nodes
                nodes = dict( [ ( n, { n } ) for n in graph.nodes_iter( ) ] )
                # collect reaction information
                reactions = defaultdict( lambda: { 0 : { "reactant" : set( ), "product" : set( ) }, 1 : set( ), -1 : set( ) } )
                recols = defaultdict( lambda: defaultdict( int ) )
                polcoll = defaultdict( set )
                newnode = ""
                for n1, n2, att in graph.edges_iter( data = True ):
                    if att[ "polarity" ] == 0:
                        reactions[ att[ "r_id" ] ][ 0 ][ "reactant" ].add( n1 )
                        reactions[ att[ "r_id" ] ][ 0 ][ "product" ].add( n2 )
                        # check for homodimer formation or dissociation
                        multimerization = False
                        if att[ "type" ] in ( "STATE_TRANSITION", "COMPLEX_ASSOCIATION" ) and graph.node[ n1 ][ "homodimer" ] < graph.node[ n2 ][ "homodimer" ]:
                            mono, mult = n1, n2
                            multimerization = True
                        elif att[ "type" ] in ( "STATE_TRANSITION", "DISSOCIATION" ) and graph.node[ n1 ][ "homodimer" ] > graph.node[ n2 ][ "homodimer" ]:
                            mono, mult = n2, n1
                            multimerization = True
                        else:
                            multimerization = False
                        if multimerization:
                            mocnt = graph.node[ mono ][ "homodimer" ]
                            mucnt = graph.node[ mult ][ "homodimer" ]
                            for i in xrange( 2 * mocnt, mucnt + 1, mocnt ):  # start with double amount of monomer b/c it is no multimer otherwise
                                newnode = "{}_{}" . format( mono, i // mocnt)
                                nodes[ mono ].add( newnode )
                    else:
                        reactions[ att[ "r_id" ] ][ att[ "polarity" ] ].add( ( n1, n2  ) )

                for r, pol in reactions.iteritems( ):
                    if len( pol[ -1 ] ) + len( pol[ 1 ] ) > 1:
                        rlabel = r + "a"
                    else:
                        rlabel = r
                    processed = False
                    for sign in ( -1, 1 ):
                        reac = dict( [ ( n, -sign ) for n in reactions[ r ][ 0 ][ "reactant" ] ] )
                        prod = dict( [ ( n, sign ) for n in reactions[ r ][ 0 ][ "product" ] ] )
                        for n1, n2 in pol[ sign ]:  # check modulators first
                            upd = dict( reac.items( ) + prod.items( ) + [ ( n1, -sign ), ( n2, sign ) ] )  # stuff to the right overwrites previous values
                            recols[ rlabel ].update( upd )
                            rlabel = rlabel[ :-1 ] + chr( ord( rlabel[ -1 ] ) + 1 )  # this one will have no consequences if only one reaction is added
                            # If they are identical, the value for n2 will have prevailed (s.a.) -> do not flag the node for the first value.
                            polcoll[ -sign ] |= set( reac ) - set( prod ) | { n1 } - { n2 }
                            polcoll[ sign ] |= set( prod ) | { n2 }
                            processed = sign == 1 # mark processed only for positive modulators
                    if not processed:  # no positive modulators means the reaction must be included once independent of modulators
                        recols[ rlabel ].update( reac.items( ) + prod.items( ) )
                        polcoll[ -1 ] |= set( reac ) - set( prod )  # same reasoning as above
                        polcoll[ 1 ] |= set( prod )

                # introduce accessory import and export reactions (to meet FBA requirements)
                templ = { 1 : "im_{}", -1 : "ex_{}" }
                for sign in ( -1, 1 ):
                    for n in set( nodes ) - polcoll[ sign ]:
                        recols[ templ[ sign ].format( n ) ][ n ] = sign

                recolsort = sorted( recols )#, key = lambda x: min( [ int( r.strip( ascii_letters ) ) for r in x.split( "+" ) ] ) )
                stream.write( "\t" . join( [ "0", "0" ] + [ str( i ) for i in xrange( 1, len( recolsort ) + 1 ) ] ) + "\n" )
                stream.write( "\t" . join( [ "CellDesigner", "Cytoscape" ] + recolsort ) + "\n" )
                for n in sorted( nodes ):
                    elements = [ graph.node[ n ][ "name" ], n ] + [ str( recols[ re ][ n ] ) for re in recolsort ]
                    for rn in sorted( nodes[ n ] ):
                        elements[ 1 ] = rn
                        stream.write( "\t" . join( elements ) + "\n" )
                stream.close( )


############################################################
#### NOT SO PUBLIC                                      ####
############################################################

    def _reset( self ):
        """reset internal storages"""
        self.graph = networkx.MultiDiGraph( shells = 0 )
        self._symbol_mapper = deepcopy( NetworkCreator.__symbol_mapper )


    def _resort_mapper( self ):
        """resort all entries in the internal mapper dictionary"""
        for s in self._symbol_mapper[ "node" ]:
            for k, l in self._symbol_mapper[ "node" ][ s ].iteritems( ):
                self._symbol_mapper[ "node" ][ s ][ k ] = sorted( set( l ) )


    def _add( self, node_def, g = None, shell = None, seeds = False, gfilter = None, complex_filter = "all" ):
        """unalias and annotate input, then add nodes"""

        if not node_def:
            return dict( )
        elif type( node_def ) != set:
            self._alert( "Wrong type for parameter(s) {:s}, {!r} expected." . format( " and " . join( [ "1" ] ), set.__name__ ) )
            return False
        if g == None:
            g = self.graph
        if shell == None:
            shell = 0
        else:
            try:
                shell = int( shell )
            except ValueError:
                #return False
                shell = 0
        if type( gfilter ) == set or gfilter == None:
            pass
        elif type( gfilter ) in ( tuple, list, dict ):
            gfilter = set( gfilter )
        elif type( gfilter ) == str:
            gfilter = { gfilter }
        else:
            self._alert( "Wrong type for parameter {!r}, {!r} expected." . format( "gfilter", "iterable" ) )
            return False

        # check and introduce pseudoreaction counter, if necessary
        if "psre" not in g.graph:
            g.graph[ "psre" ] = 0

        db_ids = zip( *node_def )[ 0 ]

        # we filter out nodes that are not wanted here, which means not adding them to the to graph
        if gfilter:
            # decide how complexes will be checked: does one matching subunit merit inclusion, or must all subunits match?
            if complex_filter in ( "all", True ):
                _decide = lambda s, gf: len( set( s ) ) == len( set( s ) & gf )  # must match all subunits
            else:
                _decide = lambda s, gf: set( s ) & gf  # must match any subunit
            annot = dict( [ ( k, v ) for k, v in self._db.query_for_annotations( db_ids ).iteritems( ) if _decide( v[ "symbol" ].split( ":" ), gfilter ) ] )
        else:
            annot = self._db.query_for_annotations( db_ids )

        grpn = dict( x = 0, y = 0, fill = '#000000', outline_width = 1.0 )  # default node graphics
        e = { "type" : "COMPLEX_ASSOCIATION", "polarity" : 0, "graphics" : dict(  target_arrow = 9, type = "line", width = 4, fill = '#880088' ) }
        s = self.getSpeciesRestriction( )  # second filter - does the species match the set restriction?
        for k, t in node_def:
            try:
                symbol = annot[ k ][ "symbol" ]
            except KeyError:
                continue
            else:
                tag = symbol
                if s and annot[ k ][ "species" ] not in s:
                    del annot[ k ]
                    continue
                d = dict( shell = shell, frontier = True, seed = int( seeds ), homodimer = 1, graphics = grpn,
                            compartment = "nucleus" )  # TODO: quick'n'dirty hack to accomodate TRANSFAC expansion
                acc_d = deepcopy( d )
                acc_d[ "class" ] = "PROTEIN"
                if t == "mirna" or ( NetworkCreator.__mirna_regexp.match( symbol.lower( ) ) and not NetworkCreator.__hostgene_regexp.match( symbol.lower( ) ) ):
                    tag += "#mir"
                    d[ "class" ] = "ANTISENSE_RNA"
                    d.update( hgnc_symbol = symbol )
                elif t == "complex":
                    tag += "#com"
                    d[ "class" ] = "COMPLEX"
                    components = sorted( symbol.split( ":" ) )
                    d.update( members = "|" . join( components ), subunits = ":" . join( components ) )
                    if len( set( components ) ) < len( components ) and len( set( components ) ) == 1:
                        d[ "homodimer" ] = len( components )
                else:
                    tag += "#pro"
                    d[ "class" ] = "PROTEIN"
                    d.update( hgnc_symbol = symbol )
                g.add_node( tag, label = tag, tag = tag, name = symbol, **d )
                self._symbol_mapper[ "node" ][ symbol.lower( ) ][ "base" ].append( tag )

                # add supplementary complex subunits (and track reactions with pseudoreaction counter)
                if t == "complex":
                    g.graph[ "psre" ] += 1
                    for c in set( components ):
                        try:
                            acc = self._symbol_mapper[ "node" ][ c.lower( ) ][ "base" ][ 0 ]
                        except IndexError:
                            acc = c  + "#pro"
                            acc_d.update( label = acc, tag = acc, hgnc_symbol = c, name = c )
                            g.add_node( acc, **acc_d )
                            self._symbol_mapper[ "node" ][ c.lower( ) ][ "base" ].append( acc )
                        e.update( instance = "{} -- {}" . format( acc, tag ), r_id = "psre{}" . format( g.graph[ "psre" ] ) )
                        g.add_edge( acc, tag, e[ "instance" ], **e )

        self._resort_mapper( )

        return annot


    def _connect( self, restrict = { "fish" : None }, final = { }, add_shells = 1, g = None, gfilter = None, cfilter = None ):
        """link nodes, optionally expanding the network with the given type(s) of molecule, working only on frontier nodes"""

        if add_shells < 0:
            return None

        if g == None:
            g = self.graph
        elif g.__class__ != networkx.MultiDiGraph:
            self._extalert( "Wrong type of argument at position {:d}: expected 'MultiDiGraph' instead of {!r}." . format( 1, g.__class__.__name__ ) )
            return False

        # check restrict skeletons
        if "fish" not in restrict:
            restrict[ "fish" ] = None
        if not final:
            final = restrict
        if "fish" not in final:
            final[ "fish" ] = None

        # prepare graph
        if "shells" not in g.graph:
            if not g.node:
                g.graph[ "shells" ] = 0
            else:
                g.graph[ "shells" ] = 1
        for n in g.nodes_iter( ):
            if "frontier" not in g.node[ n ]:
                g.node[ n ][ "frontier" ] = True
            if "shell" not in g.node[ n ]:
                g.node[ n ][ "shell" ] = g.graph[ "shells" ]

        annot = { }
        node_set = set( )
        edge_set = set( )
        grpe_0 = dict( target_arrow = 9, type = "line", width = 4, fill = '#00cccc' )
        grpe_p = dict( target_arrow = 12, type = "line", width = 4, fill = '#00cccc' )
        grpe_n = dict( target_arrow = 15, type = "line", width = 4, fill = '#00cccc' )
        while add_shells > -1:
            # The outer_shell computation assumes that all complex subunits are also available as monomeric nodes (that have an hgnc_symbol attribute).
            outer_shell = dict( [ ( n, a[ "hgnc_symbol" ] ) for n, a in g.nodes_iter( data = True ) if a[ "frontier" ] and "hgnc_symbol" in a ] )
            ref = self._db.query_for_references( outer_shell.values( ), { "alias_types" : [ "hgnc.symbol" ] }, invert = True )[ 0 ]
            ref.update( unmappable_identifier = int( 1e9 ) )  # introduce unmappable identifier for complex subunits
            # The condition in next comprehension excludes non-annotated nodes in mixed species sources from raising an exception.
            ref = dict( [  ( ref[ v.lower( ) ], k ) for k, v in outer_shell.iteritems( ) if v.lower( ) in ref  ] )

            if add_shells == 0:  # in the last round, ...
                restrict = final  # ... stop adding new nodes, but loop once more to interconnect the outermost layer
                gfilter = outer_shell.values( )
                cfilter = "all"  # no new complexes unless all subunits are already present
            else:
                g.graph[ "shells" ] += 1

            new_nodes = set( )
            new_edges = dict( )
            res = self._db.query_for_interactions( ref, restrict = restrict )
            node_set |= set( ref )
            for i in xrange( len( res ) ):
                k1, k2, t1, t2 = res[ i ][ :4 ]
                edge = k1, k2, res[ i ][ 4 ][ "database" ], res[ i ][ 4 ][ "PMIDs" ]
                if { k1, k2 } <= node_set:
                    if edge not in edge_set:
                        new_edges[ edge ] = i
                    continue
                elif restrict[ "fish" ] == set( ):  # no new node allowed -> no new connection
                    continue
                elif k1 in node_set:
                    missing = ( k2, t2 )
                else:
                    missing = ( k1, t1 )

                m = missing
                if m[ 0 ] > 1e9:
                    m = m[ 0 ], "complex"
                new_nodes.add( m )
                if edge not in edge_set:
                    new_edges[ edge ] = i


            annot.update( self._add( new_nodes, shell = g.graph[ "shells" ], gfilter = gfilter, complex_filter = cfilter ) )
            if add_shells > 0 and not annot:
                self._alert( "Something went wrong while adding new nodes." )
                return -1

            # make connections
            for i in new_edges.values( ):
                r = res[ i ]
                for j in ( 0, 1 ):
                    try:
                        symbol = annot[ r[ j ] ][ "symbol" ].lower( )
                    except KeyError:
                        try:
                            r[ j ] = ref[ r[ j ] ]
                        except KeyError:  # consistent with not adding the node in _add, either due to gene or species filtering
                            r[ j ] = False
                    else:
                        try:
                            r[ j ] = self._symbol_mapper[ "node" ][ symbol ][ "base" ][ 0 ]
                        except IndexError:
                            r[ j ] = chain( *self._symbol_mapper[ "node" ][ symbol ].values( ) ).next( )
                if False in r:
                    continue
                g.graph[ "psre" ] += 1
                e = dict( type = "DATABASE", r_id = "psre{}" . format( g.graph[ "psre" ] ) )
                e.update( r[ 4 ] )
                if e[ "database" ].lower( ) in ( "mirtarbase", "mirtarbase6", "tarbase", "mirecords" ):
                    e.update( instance = "{} -I {}" . format( r[ 0 ], r[ 1 ] ), polarity = -1, graphics = grpe_n )
                elif e[ "database" ].lower( ) in ( "regphos", ):
                    e.update( instance = "{} -o {}" . format( r[ 0 ], r[ 1 ] ), polarity = 1, graphics = grpe_p )
                else:
                    e.update( instance = "{} -D {}" . format( r[ 0 ], r[ 1 ] ), polarity = 0, graphics = grpe_0 )
                g.add_edge( r[ 0 ], r[ 1 ], e[ "instance" ], **e )
            edge_set |= set( new_edges )

            # unset frontier property
            for key in outer_shell:
                g.node[ key ][ "frontier" ] = False

            if not new_nodes:
                add_shells = -1  # no new nodes -> stop iterating
            else:
                add_shells -= 1  # next shell

        return


    def _read_CSV( self, path ):
        """convert an existing network from a CSV file into our representation"""

        return  # NOTE: not tested in this version

        try:
            stream = open( path )
        except IOError:
            self._extalert( "Error reading {!r}." . format( path ) )

        nodes = dict( )
        linecounter = 0

        for line in stream:
            linecounter += 1
            try:
                ( source, inter, target ) = line.split( )
            except ValueError:
                self._extalert( "Skipping incompatible line {:d} in {!r}." . format( linecounter, path ) )
                continue
            try:
                nodes[ source ].addInteraction( target, "molecule", inter, None )
            except KeyError:
                nodes[ source ] = NetworkNode( source, "molecule" )
                nodes[ source ].addInteraction( target, "molecule", inter, None )

        stream.close( )
        return nodes


    def _read_SBML( self, path ):
        pass


    def _read_XGMML( self, file_list, nodes ):
        """build network representation from an XGMML file using lxml's etree """

        return  # NOTE: not tested in this version

        for f in file_list:
            try:
                g = etree.parse( f )
            except ( IOError, lxml.etree.XMLSyntaxError ):
                continue

            for el in g.getroot( ):
                if el.tag[ -4 : ] == "node":
                    try:
                        ( name, molType ) = el.attrib[ "id", "type" ]
                    except KeyError:
                        continue
                    if name not in nodes.keys( ):
                        nodes[ name ] = NetworkNode( name, molType )

                elif el.tag[ -4: ] == "edge":
                    pass



