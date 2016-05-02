
#from libsbml import *
from __future__ import print_function
#from libsbml import SBMLDocument, SBMLWriter, LIBSBML_OPERATION_SUCCESS
from miRNexpander.mWBBase import mWBBaseClass
import networkx as nx
import matplotlib.pyplot as plt
#import SBOTerms
import os, sys
import time
import re

#sys.dont_write_bytecode = True
#''' The necessity for the above line suddenly popped up after libsbml had been compiled with '--enable-layout'.
#    Some bug in either python or libsbml leads to a glibc "double free or corruption" error once the python interpreter ( or debugger ) exits.
#    The script itself is executed flawlessly, however. Must be magic.
#    CAVE: Not sure yet whether dont_write_bytecode truly does the trick.'''

class NetworkWriter( mWBBaseClass ):
    """represent a network and provide the tools for transforming it into an SBML file"""

    def __init__( self, graph ):
        self.Graph = graph
        self.Nodes = { }  # TODO: remove

        try:
            self.ModelName = self.Graph.graph[ 'name' ]
        except KeyError:
            self.ModelName = "New_Network"
            self._alert( "Invalid model name - replaced with {!r}." . format( self.ModelName ) )

        self.sbmldoc = SBMLDocument( 3, 1 )
        #sbmlns = SBMLNamespaces( 3, 1, "layout", 1 )  # include graphic definitions by referring to "layout" -> Cytoscape does not evaluate them, however.
        #self.sbmldoc = SBMLDocument( sbmlns )
        #self.sbmldoc.setPkgRequired( "render", False )
        self._output_dir = ""

        self._id_replace = re.compile( "[^A-Za-z0-9_]" )  # defines characters that are not allowed in SBML identifiers

        self._output_formats = {
                "XGMML": {
                    "extension" : "xgmml",  # the file extension
                    "mirna_height" : 35,
                    "mirna_width" : 95,
            		"mirna_fill" : "#6666ff",  # hex rgb value of the node fill color for miRNAs
            		"mirna_label" : "#000000",  # hex rgb value of the node fill color for miRNAs
                    "mirna_shape" : "PARALLELOGRAM",
                    "protein_height" : 35,
                    "protein_width" : 55,
            		"protein_fill" : "#ff6660",  # hex rgb value of the node fill color for proteins
            		"protein_label" : "#000000",  # hex rgb value of the node fill color for proteins
                    "protein_shape" : "ELLIPSE",
                    "seed_height" : 35,
                    "seed_width" : 95,
            		"seed_fill" : "#000000",  # hex rgb value of the node fill color for seed
            		"seed_label" : "#00ffff",  # hex rgb value of the node label color for seeds
                    "seed_shape" : "DIAMOND",
					"activ_color" : "#00ff66",
					"activ_source" : 0,
					"activ_target" : 6,
					"activ_width" : 2,
					"inhib_color" : "#ff0066",
					"inhib_source" : 0,
					"inhib_target" : 15,
					"inhib_width" : 1,
					"symm_color" : "#444444",
					"symm_source" : 12,
					"symm_target" : 12,
					"symm_width" : 2,
					"control_color" : "#999999",
					"control_source" : 0,
					"control_target" : 12,
					"control_width" : 2,
                    },
                }


    def setOutputDir( self, target ):
        """store (and create, if not existing) the path where files will be put"""
        if target == "":
            target = "."  # fail-safe measure
        if not os.path.isdir( target ):
            try:
                os.mkdir( target )
            except OSError as e:
                self._extalert( "Unable to create directory '{}'." . format( target ) )
                return 1
        self._output_dir = target


    def buildMolID( self, molType, molName ):
        return "{}_{}" . format( molType.upper( ), re.sub( self._id_replace, "_", molName ) )


    def saveImage( self, nrepr = "spring", ft = "png", filename = "", graph = None ):
        """store an image representation of the network"""

        if graph == None:
            if self.Graph == None:
                self._alert( "No network given, not drawing." )
                return False
            else:
                graph = self.Graph

        if self._output_dir == "":
            self._extalert( "Target directory needs to be assigned with {!r} first!" )
            return 1
        if ft not in ( "png" ):
            ft = "png"
        if filename == "":
            filename = self.ModelName + "." + ft
        filename = os.path.join( self._output_dir, filename )

        if nrepr == "shells":
            shells = [ [ key for key in graph.nodes_iter( ) if graph.node[ key ][ "shell" ] == s ] for s in xrange( graph[ "shell" ] + 1 ) ]
            pos = nx.shell_layout( graph, shells )
        else:
            pos = nx.spring_layout( graph )

        plt.figure( figsize = ( 8, 8 ) )
        nx.draw( graph, pos )
        plt.savefig( filename )


    def createSBML( self, filename = "" ):
        """store a network representation in SBML format"""

        g = self.Graph
        if not g.number_of_nodes( ):
            self._alert( "No nodes were initialized, not building." )
            return 1
        if self._output_dir == "":
            self._extalert( "Target directory needs to be assigned with {!r} first!" . format( "setOutputDir" ) )
            return 1
        if filename == "":
            filename = self.ModelName + ".xml"
        filename = os.path.join( self._output_dir, filename )

        model = self.sbmldoc.createModel( )
        model.setId( self.ModelName )

        #graphics = (  model.getPlugin( "layout" )  )
        #layout = graphics.createLayout( )

        comp = model.createCompartment( )
        comp.setName( "cytosol" )
        comp.setId( "comp_CYTOSOL" )
        comp.setSize( 1e-14 )  # in liters = 1e-3 cubic meters; this is negligible for our purposes
        comp.setConstant( True )

        node_prefix = "n"
        for n in g.nodes_iter( ):
            #self._spill( "Processing node {} for inclusion in SBML file." . format( node.getName( ) ) )

            molID = node_prefix + str( n )
            node = g.node[ n ]
            sp = model.createSpecies( )
            sp.setId( molID )
            sp.setName( node[ "symbol" ] )
            sp.setCompartment( comp.getId( ) )
            sp.setHasOnlySubstanceUnits( True )
            sp.setBoundaryCondition( True )
            sp.setConstant( True )

        i = 0
        for edge in g.edges_iter( data = True ):
            i += 1  # increase reaction id counter'''
            s = node_prefix + str( edge[ 0 ] )
            t = node_prefix + str( edge[ 1 ] )
            rxn = model.createReaction( )
            rxn.setId( 'rxn_' + `i` )
            rxn.setName( edge[ 2 ][ "source_alias" ] + '->' + edge[ 2 ][ "target_alias" ] )
            rxn.setFast( False )
            rxn.setReversible( False )
            reactant = rxn.createReactant( )
            reactant.setSpecies( s )
            reactant.setConstant( True )
            product = rxn.createProduct( )
            product.setSpecies( t )
            product.setConstant( False )

        modelerrors = self.sbmldoc.checkInternalConsistency( )
        if modelerrors:
            self._extalert( "SBML inconsistencies detected." )
            self._spill( self.sbmldoc.printErrors( ) )

        writer = SBMLWriter( )
        try:
            writer.writeSBMLToFile( self.sbmldoc, filename )
            self._spill( "Finished writing SBML file {!r}." . format( filename ) )
        except:
            self._extalert( "Cannot write SBML to file {!r}." . format( filename ) )


    def createCytoscape( self, filename = "" ):
        """store a network representation in Cytoscape-readable format"""

        if not self.Nodes or type( self.Nodes ) is not type( { } ):
            self._alert( "{}: No nodes were initialized, not building." . format( "createCytoscape( )" ) )
            return 1
        if self._output_dir == "":
            self._extalert( "Target directory needs to be assigned with {!r} first!" . format( "setOutputDir" ) )
            return 1
        if filename == "":
            filename = self.ModelName + ".csv"
        filename = os.path.join( self._output_dir, filename )

        output = [ ]
        for node in self.Nodes.values( ):
            #print( "Processing node {} for inclusion in Cytoscape file." . format( node.getName( ) ) )

            source_type = node.getMolType( )
            for link in node.getLinks( ):
                target_type = link.getTarget( ).getMolType( )
                reg = link.getEffect( )
                # edit the following line to change the column separator (tabstop is standard)
                if link.getIsDirected( ):
                    interaction_string = " -> " . join( [ source_type, reg, target_type ] )
                else:
                    interaction_string = source_type + " | " + reg + " | " + target_type
                output.append( "{}\t{}\t{}\n" . format( node.getName( ), interaction_string, link.getTarget( ).getName( ) ) )

        try:
            stream = open( filename, "w" )
        except IOError as e:
            self._extalert( "{}" . format( e.strerror ) )
            return 1
        stream.write( "" . join( output ) )
        stream.close( )
        self._spill( "Finished writing Cytoscape file {!r}." . format( filename ) )


    def createGML( self, filename = "", graph = None, gml_dict = None ):
        """store a network representation in GML format"""

        if graph == None:
            if self.Graph == None:
                self._alert( "No network given, not drawing." )
                return False
            else:
                g = self.Graph.copy( )  # copying here is important because we will make changes to the graph solely for output purposes

        if not gml_dict:
            gml_dict = self._output_formats[ "XGMML" ]  # re-use xgmml for now
        if self._output_dir == "":
            self._extalert( "Target directory needs to be assigned with {!r} first!" . format( "setOutputDir" ) )
            return False
        if filename == "":
            filename = self.ModelName + '.gml'
        filename = os.path.join( self._output_dir, filename )

        # before writing, we relabel the nodes with a meaningful name and add graphic stuff
        labels = { }
        for n in g.nodes_iter( ):
            labels[ n ] = g.node[ n ][ "symbol" ]
            if g.node[ n ][ "seed" ]:
                t = "seed"
            else:
                t = g.node[ n ][ "molType" ]
            try:
                g.node[ n ][ "graphics" ] = {
                    "w" : gml_dict[ t + "_width" ],
                    "h" : gml_dict[ t + "_height" ],
                    "fill" : gml_dict[ t + "_fill" ],
                    "type" : gml_dict[ t + "_shape" ],
                    #"outline"
                    #"outline_width"
                    }
            except KeyError:
                self._alert( "Could not find graphical parameters for node {!r} of type {!r}." . format( n, t ) )
        g = nx.relabel_nodes( g, labels )

        # do similar stuff for edges
        for e in g.edges_iter( data = True ):
            reg = e[ 2 ][ "regulation" ]
            if reg in ( "inhibition", "repression", "suppression", "downregulation", "intranslation", "sequestration" ):
                t = "inhib"
            elif reg in ( "activation", "expression", "derepression", "upregulation", "translation", "liberation" ):
                t = "activ"
            elif reg in ( "regulation", "control", "derivative", "inclusion" ):  # if the relation is unknown
                t = "control"
            else:  # symmetric or undefined relationship, e.g. between two protein binding partners
                t = "symm"
            e[ 2 ][ "graphics" ] = {
                    "type" : "line",
                    "fill" : gml_dict[ t + "_color" ],
                    #"Line" : {  "w" : gml_dict[ t + "_width" ]  },
                    #"Line" : { },
                    "source_arrow" : gml_dict[ t + "_source" ],
                    "target_arrow" : gml_dict[ t + "_target" ],
                    }

        nx.write_gml( g, filename )
        self._spill( "Finished writing GML file {!r}." . format( filename ) )


    def createXGMML( self, filename = "", xgmml_dict = None ):
        """store a network representation in XGMML format"""

        ### internal subroutine ###
        def _convert_node( node, gs ):
            """return an XGMML node string"""

            nodeID = node[ 0 ]  # might change in the future
            node = node[ 1 ]
            node_markup = '{i1}<node label="{0}" id="{1}">\n'
            #for k in set( node ) - set( [ "id" ] ):
            for k in set( node ):
                if type( node[ k ] ) == str:
                    t = "string"
                    v = node[ k ]
                elif type( node[ k ] ) == int:
                    t = "integer"
                    v = str( node[ k ] )
                elif type( node[ k ] ) == bool:
                    t = "boolean"
                    v = str( int( node[ k ] ) )
                node_markup += '    <att name="{}" type="{}" value="{}"/>\n' . format( k, t, v )
            node_markup += ( '{i2}<graphics type="{gnt}" fill="{gnf}" h="{gnh}" w="{gnw}" width="{gnow}" outline="{gnoc}" '
                            + 'cy:nodeLabel="{0}" cy:nodeLabelColor="{glc}"/>\n'
                            + '{i1}</node>\n'
                            )
            nodeName = node[ "symbol" ]
            #nodeID = node[ "id" ]  # might change in the future
            nodeType = node[ "molType" ]
            if node[ "seed" ]:
                node_emph = gs[ "seed_graphics" ]
            else:
                try:
                    node_emph = gs[ nodeType.lower( ) + "_graphics" ]
                except:
                    if nodeType == "pri-miRNA":
                        node_emph = gs[ "pri-mirna_graphics" ]
                    elif nodeType == "miRNA":
                        node_emph = gs[ "mirna_graphics" ]
                    else:
                        node_emph = gs[ "prot_graphics" ]

            format_dict = dict( gs[ "base_values" ].items( ) + gs[ "node_graphics" ].items( ) + node_emph.items( ) )
            #node_string = gs[ "node_markup" ] . format( nodeName, nodeID, "#n/a#", nodeType, **format_dict )
            node_string = node_markup . format( nodeName, nodeID, **format_dict )
            return node_string


        ### internal subroutine ###
        def _convert_edge( edge, mapper, edge_counter, gs ):
            """return an XGMML edge string"""

            edge_markup = '{i1}<edge label="{0}" source="{1}" target="{2}" id="{3}" cy:directed="1">\n'
            #for k in set( edge[ 2 ] ) - set( [ "id", "source", "target" ] ):
            for k in set( edge[ 2 ] ) - set( [ "source", "target" ] ):
                edge_markup += '    <att name="{}" type="string" value="{}"/>\n' . format( k, edge[ 2 ][ k ] )
            edge_markup += ( '{i2}<graphics fill="{gec}" width="{gew}" cy:sourceArrow="{gsa}" cy:targetArrow="{gta}"/>\n'
                            + '{i1}</edge>\n'
                            )

            # map non-numerical ids to changed ones
            try:
                source = change[ edge[ 0 ] ]
            except KeyError:
                source = edge[ 0 ]
            try:
                target = change[ edge[ 1 ] ]
            except KeyError:
                target = edge[ 1 ]

            label = " -&gt; "
            reg = edge[ 2 ][ "regulation" ]
            if reg in ( "inhibition", "repression", "suppression", "downregulation", "intranslation", "sequestration" ):
                edge_emph = gs[ "inhib_graphics" ]
            elif reg in ( "activation", "expression", "derepression", "upregulation", "translation", "liberation" ):
                edge_emph = gs[ "activ_graphics" ]
            elif reg in ( "regulation", "control", "derivative", "inclusion" ):  # if the relation is unknown
                edge_emph = gs[ "control_graphics" ]
            else:  # symmetric or undefined relationship, e.g. between two protein binding partners
                edge_emph = gs[ "symm_graphics" ]
                label = " &lt;-&gt; "
            label = edge[ 2 ][ "source_alias" ] + label + edge[ 2 ][ "target_alias" ]

            format_dict = dict( gs[ "base_values" ].items( ) + edge_emph.items( ) )
            edge_string = edge_markup . format( label, source, target, edge_counter, **format_dict )
            return edge_string


        ### code proper of 'createXGMML( )' starts here
        if not self.Graph:
            self._extalert( "{}: No graph initialized, not building." )
            return False
        else:
            g = self.Graph
        if self._output_dir == "":
            self._extalert( "Target directory needs to be assigned with {!r} first!" . format( "setOutputDir" ) )
            return False
        if not xgmml_dict:
            xgmml_dict = self._output_formats[ "XGMML" ]
        if filename == "":
            filename = self.ModelName + '.' + xgmml_dict[ "extension" ]
        filename = os.path.join( self._output_dir, filename )

        # constant definitions
        indent_str = "  "
        max_indent = 4
        try:
            graph_settings = {  # this dictionary stores values that need to be available to subroutines
                "node_graphics" : {
                        "gnow" : "1",  # not sure yet what that means; outline thickness ?
                        "gnoc" : "#c0c0c0",  # outline color (?)
                        "gns"  : 0,  # was this node a seed for the network?
                        },
                "pri-mirna_graphics" : {
                        "gnh" : "{:.1f}" . format( xgmml_dict[ "pri-mirna_height" ] ),  # node height
                        "gnw" : "{:.1f}" . format( xgmml_dict[ "pri-mirna_width" ] ),  # node width
                        "gnt" : xgmml_dict[ "pri-mirna_shape" ],  # node shape
                        "gnf" : xgmml_dict[ "pri-mirna_fill" ],  # node color
                        "glc" : xgmml_dict[ "pri-mirna_label" ],  # label color
                        },
                "mirna_graphics" : {
                        "gnh" : "{:.1f}" . format( xgmml_dict[ "mirna_height" ] ),  # node height
                        "gnw" : "{:.1f}" . format( xgmml_dict[ "mirna_width" ] ),  # node width
                        "gnt" : xgmml_dict[ "mirna_shape" ],  # node shape
                        "gnf" : xgmml_dict[ "mirna_fill" ],  # node color
                        "glc" : xgmml_dict[ "mirna_label" ],  # label color
                        },
                "protein_graphics" : {
                        "gnh" : "{:.1f}" . format( xgmml_dict[ "protein_height" ] ),  # node height
                        "gnw" : "{:.1f}" . format( xgmml_dict[ "protein_width" ] ),  # node width
                        "gnt" : xgmml_dict[ "protein_shape" ],  # node shape
                        "gnf" : xgmml_dict[ "protein_fill" ],  # node color
                        "glc" : xgmml_dict[ "protein_label" ],  # label color
                        },
                "seed_graphics" : {
                        "gnh" : "{:.1f}" . format( xgmml_dict[ "seed_height" ] ),  # node height
                        "gnw" : "{:.1f}" . format( xgmml_dict[ "seed_width" ] ),  # node width
                        "gnt" : xgmml_dict[ "seed_shape" ],  # node shape
                        "gnf" : xgmml_dict[ "seed_fill" ],  # node color
                        "glc" : xgmml_dict[ "seed_label" ],  # label color
                        "gns" : 1,
                        },
                "activ_graphics" : {
                        "gec" : xgmml_dict[ "activ_color" ],  # edge color
                        "gew" : "{:.1f}" . format( xgmml_dict[ "activ_width" ] ),  # edge thickness
                        "gsa" : "{:.0f}" . format( xgmml_dict[ "activ_source" ] ),  # arrow end at source
                        "gta" : "{:.0f}" . format( xgmml_dict[ "activ_target" ] ),  # arrow end at target
                        },
                "inhib_graphics" : {
                        "gec" : xgmml_dict[ "inhib_color" ],  # edge color
                        "gew" : "{:.1f}" . format( xgmml_dict[ "inhib_width" ] ),  # edge thickness
                        "gsa" : "{:.0f}" . format( xgmml_dict[ "inhib_source" ] ),  # arrow end at source
                        "gta" : "{:.0f}" . format( xgmml_dict[ "inhib_target" ] ),  # arrow end at target
                        },
                "control_graphics" : {
                        "gec" : xgmml_dict[ "control_color" ],  # edge color
                        "gew" : "{:.1f}" . format( xgmml_dict[ "control_width" ] ),  # edge thickness
                        "gsa" : "{:.0f}" . format( xgmml_dict[ "control_source" ] ),  # arrow end at source
                        "gta" : "{:.0f}" . format( xgmml_dict[ "control_target" ] ),  # arrow end at target
                        },
                "symm_graphics" : {
                        "gec" : xgmml_dict[ "symm_color" ],  # edge color
                        "gew" : "{:.1f}" . format( xgmml_dict[ "symm_width" ] ),  # edge thickness
                        "gsa" : "{:.0f}" . format( xgmml_dict[ "symm_source" ] ),  # arrow end at source
                        "gta" : "{:.0f}" . format( xgmml_dict[ "symm_target" ] ),  # arrow end at target
                        },
                "base_values" : dict(
                            [ ( "nn", self.ModelName ), ( "date", time.strftime( "%Y-%m-%d %H:%M:%S", time.localtime( ) ) ) ]
                            + [ ( "i" + `i`, indent_str * i ) for i in xrange( 1, max_indent + 1 ) ]
                        )
            }
        except KeyError as e:
            self._extalert( "Error in the XGMML format dictionary.\n{}" . format( e ) )
            return 1

        previousNodes = { }
        output = [ str( '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
            '<graph label="{nn}" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" '
                'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:cy="http://www.cytoscape.org" '
                'xmlns="http://www.cs.rpi.edu/XGMML" directed="1">\n'
            '{i1}<att name="documentVersion" value="1.1"/>\n'
            '{i1}<att name="networkMetadata">\n'
            '{i2}<rdf:RDF>\n'
            '{i3}<rdf:Description rdf:about="http://www.cytoscape.org/">\n'
            '{i4}<dc:type>[miRNA,Protein]-Protein Interaction</dc:type>\n'
            '{i4}<dc:description>N/A</dc:description>\n'
            '{i4}<dc:identifier>N/A</dc:identifier>\n'
            '{i4}<dc:date>{date}</dc:date>\n'
            '{i4}<dc:title>{nn}</dc:title>\n'
            '{i4}<dc:source>http://www.cytoscape.org/</dc:source>\n'
            '{i4}<dc:format>Cytoscape-XGMML</dc:format>\n'
            '{i3}</rdf:Description>\n'
            '{i2}</rdf:RDF>\n'
            '{i1}</att>\n'
            '{i1}<att type="string" name="backgroundColor" value="#ffffff"/>\n' ) . format( **graph_settings[ "base_values" ] )
        ]


        # replace non-numerical ids, e.g. in miRNAs where the MIMAT is used
        fine = [ ]
        change = { }
        for node in g.nodes_iter( ):
            try:
                fine.append( int( node ) )
            except ValueError:
                change[ node ] = 0
        m = max( fine )
        for node in change:
            m += 1
            #g.node[ node ][ "id" ] = m
            change[ node ] = m

        # convert nodes and edges
        for node in g.nodes_iter( data = True ):
            #output.append( _convert_node( node[ 1 ], graph_settings ) )
            output.append( _convert_node( node, graph_settings ) )
        for edge in g.edges_iter( data = True ):
            m += 1
            output.append( _convert_edge( edge, change, m, graph_settings ) )

        output.append( "</graph>\n" )
        try:
            stream = open( filename, "w" )
        except IOError as e:
            self._extalert( "Unable to write to file '{}'." . format( filename ) )
            return 1
        stream.write( "" . join( output ) )
        stream.close( )
        self._spill( "Finished writing XGMML file {!r}." . format( filename ) )


