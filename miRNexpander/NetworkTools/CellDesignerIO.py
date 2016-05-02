#!/usr/bin/python

from __future__ import print_function

import sys  # to attach to the system
import os  # for path operations
import codecs  # required to write utf-8 output to a file
from copy import deepcopy  # for copying nested types by value
from itertools import combinations as comb, chain  # for (optionally) linking complex subunits among each other
import re  # regulary expression
from xml.dom.minidom import parse  # xml parser
import networkx  # graph tool
from numpy import median  # for calculating median of expression values
from operator import lt as oplt, gt as opgt  # for flexible indicator matrix calculation
from string import ascii_letters  # for not having to write "ABC..xyz"

from collections import defaultdict

#from miRNexpander.mWBBase import mWBBaseClass
from miRNexpander.NetworkTools.NetworkCreator import NetworkCreator


class CellDesignerIO( NetworkCreator ):

    def __init__( self, db_handler ):
        """initialization"""

        NetworkCreator.__init__( self, db_handler )

        # storage for the parsed CellDesigner file
        self.xml_object = None

        # keep track of which CellDesigner object is represented by which node in the internal MultiDiGraph
        self._representative = dict( spe = dict( ), ref = dict( ), ali = dict( ) )

        # GML node and edge attribute identifiers have a restricted character space (at least in Cytoscape)
        self._gml_invalid = re.compile( "[^A-Za-z0-9_]" )


############################################################
#### PUBLIC                                             ####
############################################################

    ### XML INPUT
    def parse_file( self, path ):
        """read a CellDesigner file into the XML object"""

        self._reset( )
        try:
            self.xml_object = parse( path )
        except IOError:
            self._spill( "Unable to read input file {!r}, aborting." . format( path ) )
            sys.exit( )
        except xml.parsers.expat.ExpatError as e:
            self._spill( "There was a parsing error, please check if {!r} is a valid XML document!" . format( path ) )
            self._spill( "Error message:\n\t" + e.str )
            sys.exit( )

    ### XML OUTPUT
    def write_pretty_xml( self, path, indent = "  " ):
        """write the internal XML object with indentation to a file"""

        try:
            self.xml_object.writexml( codecs.open( path, "w", "utf-8" ), addindent = indent )
        except IOError:
            self._spill( "Unable to open output file {!r} for writing." . format( path ) )


    def write_xml( self, path ):
        """write the internal XML object with indentation to a file"""

        try:
            self.xml_object.writexml( codecs.open( path, "w", "utf-8" ) )
        except IOError:
            self._spill( "Unable to open output file {!r} for writing." . format( path ) )


    ### CONVERSION from XML to networkx
    def xml_to_graph( self, transform = { } ):
        """create a networkx representation from the CellDesigner XML object
    - transform can have the following keys:
    -- fold_to_symbol: (arbitrary)
    -- link_subunits: (arbitrary)
    -- bridge_reactions: iterable of type identifiers for reactions to be removed, potentially skipping the input node
    -- bridge_nodes: iterable of class identifiers for nodes to be bridged and removed
"""

        # parse compartment information
        cmpt = { "" : "" }  # stores compartment names; empty default value for nodes that do not have the attribute
        for x in self.xml_object.getElementsByTagName( "compartment" ):
            i, n = [ x.getAttribute( a ) for a in ( "id", "name" ) ]
            if n:
                cmpt[ i ] = n
            else:
                cmpt[ i ] = i

        # parse species information
        sp = dict( )  # stores species
        sp_template = {
            "catalyzes" : "",
            "members" : "",
            "compartment" : "",
            "metaid" : "",
            "homodimer" : 1,
        }
        ref = dict( )  # stores references to species - links species in complexes to their annotated free counterparts

        # prepare annotation matching
        miriam = dict( self._db._sql( "SELECT `URN`, `namespace` FROM Actor_xrefs" ) )

        # generate initial species dictionary - tag "species" lists only the free (=uncomplexed) components
        for x in self.xml_object.getElementsByTagName( "species" ):

            d = {
                    "class" : x.getElementsByTagName( "celldesigner:class" )[ 0 ].childNodes[ 0 ].data,  # PROTEIN, ANTISENSE_RNA etc
                    "catalyzes" : "|" . join( [ cat.getAttribute( "reaction" ) for cat in x.getElementsByTagName( "celldesigner:catalyzed" ) ] ),
                    "members" : "",
                    "mods" : dict( [ ( mod.getAttribute( "residue" ), { "state" : mod.getAttribute( "state" ) } ) for mod in x.getElementsByTagName( "celldesigner:modification" ) ] ),
            }
            try:
                d[ "homodimer" ] = int( x.getElementsByTagName( "celldesigner:homodimer" )[ 0 ].childNodes[ 0 ].data )
            except ( IndexError, ValueError ):
                d[ "homodimer" ] = 1

            ### convert annotation
            annots = {  rdf.getAttribute( "rdf:resource" ) for rdf in x.getElementsByTagName( "rdf:li" )  }
            unknown = set( )
            annots_collected = dict( )
            # the node will have only ONE annotation per attribute, meaning only one identifier per database is kept
            for ann in [ a.split( ":" ) for a in annots ]:
                urn = ":" . join( ann[ :3 ] )
                try:
                    attrib = self._gml_invalid.sub( "_", miriam[ urn ] )
                except KeyError:
                    unknown.add( ":" . join( ann ) )
                else:
                    d.update( {  attrib : ann[ -1 ]  } )
                    try:
                        annots_collected[ attrib ].add( ann[ -1 ] )
                    except KeyError:
                        annots_collected[ attrib ] = { ann[ -1 ] }
            d.update( [ ( k + "_all", "|" . join( v ) ) for k, v in annots_collected.items( ) if len( v ) > 1 ], other_annotation = "|".  join( unknown ) )

            # assign species and molecule identifiers to each other
            d[ "species" ] = x.getAttribute( "id" )
            sref = x.getElementsByTagName( "celldesigner:{}Reference" . format( d[ "class" ].replace( "_", "" ).lower( ) ) )
            if sref:
                d[ "reference" ] = sref[ 0 ].childNodes[ 0 ].data
                ref[ d[ "species" ] ] = d[ "reference" ]
                try:
                    ref[ d[ "reference" ] ].append( d[ "species" ] )
                except KeyError:
                    ref[ d[ "reference" ] ] = [ d[ "species" ] ]
            else:
                d[ "reference" ] = ""
            d.update( [ ( a, x.getAttribute( a ) ) for a in "name", "compartment", "metaid" ] )
            d[ "compartment" ] = cmpt[ d[ "compartment" ] ]
            sp[ d[ "species" ] ] = d


        # parse additional information about species and complete species dictionary
        memspec = defaultdict( list )
        if "fold_to_symbol" not in transform:
            pass
        else:
            for x in self.xml_object.getElementsByTagName( "celldesigner:species" ):
                sid = x.getAttribute( "id" )
                # If this component is part of a complex, it needs to be newly introduced to the species dictionary.
                # It can be linked to its free species via the class reference, so we need to extract that first.
                cl = x.getElementsByTagName( "celldesigner:class" )[ 0 ].childNodes[ 0 ].data
                sref = x.getElementsByTagName( "celldesigner:{}Reference" . format( cl.replace( "_", "" ).lower( ) ) )
                try:
                    rid = sref[ 0 ].childNodes[ 0 ].data
                except IndexError:  # no reference in this node
                    rid = ""
                else:
                    ref[ sid ] = rid
                try:
                    refsid = ref[ rid ][ 0 ]  # simply take the first one, it does not matter here
                except KeyError:  # reference has no species
                    refsid = None
                if refsid:  # reference one of the appropriate species id, BUT retouch modifications!!!
                    sp[ sid ] = deepcopy( sp[ refsid ] )
                    sp[ sid ][ "mods" ] = dict( [ ( mod.getAttribute( "residue" ), { "state" : mod.getAttribute( "state" ) } ) for mod in x.                     getElementsByTagName( "celldesigner:modification" ) ] )
                    sp[ sid ][ "metaid" ] = ""  # celldesigner:species have no metaid
                    ref[ rid ].append( sid )
                else:
                    #self._spill( "No free species for id {} ({}), replacing with dummy." . format( x.getAttribute( "id" ), x.getAttribute( "name" ) ) )
                    sp[ sid ] = sp_template.copy( )
                    sp[ sid ].update( {
                        "name" : x.getAttribute( "name" ),
                        "species" : sid,
                        "class" : cl,
                        "reference" : rid,
                        "mods" : dict( [ ( mod.getAttribute( "residue" ), { "state" : mod.getAttribute( "state" ) } ) for mod in x.getElementsByTagName( "celldesigner:modification" ) ] ),
                    } )
                    try:
                        d[ "homodimer" ] = int( x.getElementsByTagName( "celldesigner:homodimer" )[ 0 ].childNodes[ 0 ].data )
                    except ( IndexError, ValueError ):
                        d[ "homodimer" ] = 1
                    if rid:
                        try:
                            ref[ rid ].append( sid )
                        except KeyError:
                            ref[ rid ] = [ sid ]

                # populate members dictionary for later use (we need to finish the iteration over all species before we can add members)
                for s in x.getElementsByTagName( "celldesigner:complexSpecies" ):
                    memspec[ s.childNodes[ 0 ].data ].append( x.getAttribute( "id" ) )

        # transform members dictionary
        for sid in memspec:
            try:
                sp[ sid ][ "members" ] = "|" . join( memspec[ sid ] )
            except KeyError:
                self._alert( "CellDesigner inconsistency: Complex {} not found in listOfSpecies." . format( sid ) )

        # read protein list to add modification positions and their labels to the corresponding species
        for x in self.xml_object.getElementsByTagName( "celldesigner:protein" ):
            d = dict( )
            for res in x.getElementsByTagName( "celldesigner:modificationResidue" ):
                d[ res.getAttribute( "id" ) ] = res.getAttribute( "name" )
            try:
                sids = ref[ x.getAttribute( "id" ) ]
            except KeyError:  # no reference to this protein among free species
                self._alert( "Protein {}'s modification residues could not be assigned to a species." . format( x.getAttribute( "id" ) ) )
                continue
            for sid in sids:
                for mid in d:
                    try:
                        sp[ sid ][ "mods" ][ mid ][ "name" ] = d[ mid ]
                    except KeyError:  # means no state for this residue was specified in the species
                        sp[ sid ][ "mods" ][ mid ] = { "name" : d[ mid ], "state" : "" }

        # create name suffix with modification information
        for s in sp.values( ):
            s[ "suffix" ] = "|" .  join( [ "_" . join( [ v[ "name" ], v[ "state" ][ :3 ] ] ) for k, v in sorted( s[ "mods" ].items( ) ) if v[ "state" ] != "" ] )
            if s[ "homodimer" ] != 1:
                if s[ "suffix" ]:
                    s[ "suffix" ] = "hm{}|{}" . format( s[ "homodimer" ], s[ "suffix" ] )
                else:
                    s[ "suffix" ] = "hm{}" . format( s[ "homodimer" ] )
            if s[ "suffix" ]:
                s[ "suffix" ] = "|" + s[ "suffix" ]

        # convert modification dict to simple string
        for sid in sp:
            if type( sp[ sid ][ "mods" ] ) == dict:  # necessary because we might iterate multiple times over one entry (because of species references)
                sp[ sid ][ "mods" ] = "|" . join( [ ":" . join( [ k, v[ "name" ], v[ "state" ] ] ) for k, v in sorted( sp[ sid ][ "mods" ].items( ) ) ] )

        # parse node (= species and complexes) positions and add nodes to the graph
        self.graph = networkx.MultiDiGraph( )
        re_complex = re.compile( "^Complex.*\(([^)]+)\)" )  # to extract complex subunits
        re_script_start = re.compile( "_su(b|per)_" )
        re_script_end = re.compile( "_endsu(b|per)_" )
        ref2 = dict( )  # will hold references to the common name of essentially identical species
        count = dict( )
        replaces = defaultdict( set )
        if "fold_to_symbol" in transform:
            dup_dict = defaultdict( lambda: dict( name = None, ref = None, sp = None ) )
            for x in self.xml_object.getElementsByTagName( "celldesigner:speciesAlias" ) + self.xml_object.getElementsByTagName( "celldesigner:complexSpeciesAlias" ):
                species = x.getAttribute( "species" )
                activity = x.getElementsByTagName( "celldesigner:activity" )[ 0 ].childNodes[ 0 ].data
                alias = x.getAttribute( "id" )
                tag = re_script_end.sub( "", re_script_start.sub( "_", sp[ species ][ "name" ].replace( "_br_", " " ).strip( ) ) )
                check_tag = ( tag.replace( " ", "" ) + "#" + sp[ species ][ "class" ].lower( )[ :3 ] ).lower( )
                if alias[ 0 ] == "c":
                    # build complex tag from species labels
                    su = sorted( [ ref2[ s ][ :-4 ] for s in memspec[ species ] ] )
                    tag = "+" . join( su )
                elif not dup_dict[ check_tag ][ "name" ]:
                    dup_dict[ check_tag ] = dict( name = sp[ species ][ "name" ], ref = sp[ species ][ "reference" ], sp = species )
                elif dup_dict[ check_tag ][ "ref" ] != sp[ species ][ "reference" ]:
                    # heuristic check for proteins that have multiple reference entries even though they are supposed to have the same
                    self._spill( "'{}' ({}) looks like '{}' ({}) but has a different reference ({} vs {})."
                                    . format( sp[ species ][ "name" ], species,
                                            dup_dict[ check_tag ][ "name" ], dup_dict[ check_tag ][ "sp" ],
                                            sp[ species ][ "reference" ], dup_dict[ check_tag ][ "ref" ] )
                    )
                tag += sp[ species ][ "suffix" ]
                if activity == "active":
                    tag += "|active"
                unique = tag + "#" + sp[ species ][ "class" ].lower( )[ :3 ]  # necessary to distinguish RNA and protein of the same name, among others
                try:
                    ref2[ alias ] = ref2[ unique ]
                except KeyError:  # means that node has not been added yet
                    # link node identifiers in the references
                    ref2[ tag ] = unique  # this is the one that can be checked early in the above process
                    ref2[ alias ] = unique  # this is the one that will be used for assigning reactions from free components
                    ref2[ species ] = unique  # this is the one that will be used for assigning reactions from complexes

                    # haul in node information and add node
                    d = {
                        "activity" : activity,
                        "graphics" : {
                                "fill" : "#" + x.getElementsByTagName( "celldesigner:paint" )[ 0 ].getAttribute( "color" )[ -6: ],
                                "outline_width" : float( x.getElementsByTagName( "celldesigner:singleLine" )[ 0 ].getAttribute( "width" ) ),
                                "x" : 0,
                                "y" : 0,
                        },
                    }
                    try:  # check if we ignore this alias due to its coloring
                        ignored_color = d[ "graphics" ][ "fill" ] in transform[ "ignore_alias_colors" ]
                    except KeyError:
                        pass
                    else:
                        if ignored_color: continue
                    d.update( sp[ species ] )
                    if alias[ 0 ] == "c":
                        d[ "alias" ] = alias
                        d[ "graphics" ][ "fill" ] = "#E0E0E0"  # make complexes a little grayer, increases visibility in Cytoscape
                    del d[ "suffix" ]
                    del d[ "species" ]
                    d[ "tag" ] = unique
                    self.graph.add_node( unique, **d )
                else:
                    ref2[ tag ] = ref2[ unique ]  # TODO: check if this is logically consistent
                    ref2[ species ] = ref2[ unique ]

                self._representative[ "spe" ][ species ] = unique
                self._representative[ "ali" ][ alias ] = unique
                try:
                    ref2[ ref[ species ] ] = unique  # this one will also be used for assigning reactions from complexes
                    self._representative[ "ref" ][ ref[ species ] ] = unique
                except KeyError:
                    pass
                replaces[ unique ].add( species )

            # add represents column to nodes
            for n, rep in replaces.iteritems( ):
                self.graph.node[ n ][ "replaces" ] = "+" . join( sorted( rep ) )

        else:
            for x in self.xml_object.getElementsByTagName( "celldesigner:speciesAlias" ) + self.xml_object.getElementsByTagName( "celldesigner:complexSpeciesAlias" ):
                d = {
                    "alias" : x.getAttribute( "id" ),
                    "species" : x.getAttribute( "species" ),
                    "activity" : x.getElementsByTagName( "celldesigner:activity" )[ 0 ].childNodes[ 0 ].data,
                    "graphics" : dict( [ ( k, float( v ) ) for k, v in x.getElementsByTagName( "celldesigner:bounds" )[ 0 ].attributes.items( ) ] ),
                }
                d[ "graphics" ].update( {
                    "fill" : "#" + x.getElementsByTagName( "celldesigner:paint" )[ 0 ].getAttribute( "color" )[ -6: ],
                    "outline_width" : float( x.getElementsByTagName( "celldesigner:singleLine" )[ 0 ].getAttribute( "width" ) ),
                } )
                try:
                    d.update( sp[ d[ "species" ] ] )
                except KeyError:
                    self._spill( "Species {} not found for node {} - looks like a CellDesigner nonconformity." . format( d[ "species" ], d[ "alias" ] ) )
                    continue
                if d[ "activity" ] == "active":
                    d[ "suffix" ] += "|active"

                # generate and prettify the displayed name
                d[ "tag" ] = re_script_end.sub( "", re_script_start.sub( "_", d[ "name" ].replace( "_br_", " " ) ) )
                if d[ "alias" ][ 0 ] == "c":
                    try:
                        d[ "tag" ] = re_complex.match( d[ "tag" ] ).group( 1 ).strip( )
                    except AttributeError:
                        pass
                    d[ "tag" ] = d[ "tag" ].replace( " / ", ":" )
                    d[ "graphics" ][ "fill" ] = "#E0E0E0"  # make complexes a little grayer, increases visibility in Cytoscape
                d[ "tag" ] += d[ "suffix" ]
                try:
                    count[ d[ "tag" ] ] += 2
                except KeyError:
                    count[ d[ "tag" ] ] = 0
                else:
                    d[ "tag" ] = d[ "tag" ].center( len( d[ "tag" ] ) + count[ d[ "tag" ] ] )

                ref2[ d[ "alias" ] ] = d[ "tag" ]
                ref2[ d[ "tag" ] ] = d[ "tag" ]
                del d[ "suffix" ]
                self.graph.add_node( d[ "tag" ], **d )
                self._representative[ "spe" ][ d[ "species" ] ] = d[ "tag" ]
                try:
                    self._representative[ "ref" ][ ref[ d[ "species" ] ] ] = d[ "tag" ]
                    # fails, e.g., with complex subunits that do not have a monomer equivalent and thus lack a reference; no reason to stop, though
                except KeyError:
                    pass
                self._representative[ "ali" ][ d[ "alias" ] ] = d[ "tag" ]

        # parse reactions
        revised_reactions = dict( )
        for x in self.xml_object.getElementsByTagName( "reaction" ):
            d = {
                    "r_id" : x.getAttribute( "id" ),
                    "metaid" : x.getAttribute( "metaid" ),
                    "type" : x.getElementsByTagName( "celldesigner:reactionType" )[ 0 ].childNodes[ 0 ].data,
                    "annotation" : "|" . join( [ ":" . join( rdf.getAttribute( "rdf:resource" ).split( ":" )[ -2: ] ) for rdf in x.getElementsByTagName( "rdf:li" ) ] ),
                    "graphics" : {
                        "type" : "line",
                        "width" : 4, #x.getElementsByTagName( "celldesigner:line" )[ 0 ].getAttribute( "width" ),
                        "fill" : "#" + x.getElementsByTagName( "celldesigner:line" )[ 0 ].getAttribute( "color" )[ -6: ],
                        "source_arrow" : 0,
                        "target_arrow" : 6,
                    },
            }
            try:  # check if we ignore this alias due to its coloring
                ignored_color = d[ "graphics" ][ "fill" ] in transform[ "ignore_reaction_colors" ]
            except KeyError:
                pass
            else:
                if ignored_color: continue
            if d[ "type" ] == "PHYSICAL_STIMULATION":  # TODO: list all possibilities
                d[ "polarity" ] = 1
            elif d[ "type" ] == "NEGATIVE_INFLUENCE":  # TODO: list all possibilities
                d[ "polarity" ] = -1
            else:
                d[ "polarity" ] = 0
            r = {
                # base reactants
                "srcs" : [ ( e.getAttribute( "alias" ), e.getAttribute( "species" ) ) for e in x.getElementsByTagName( "celldesigner:baseReactant" ) ],
                # additional reactants (same conceptual status as base reactants, but handled differently in CellDesigner)
                "adds" : [ ( e.getAttribute( "alias" ), e.getAttribute( "reactant" ) ) for e in x.getElementsByTagName( "celldesigner:reactantLink" ) ],
                # (base) products
                "tgts" : [ ( e.getAttribute( "alias" ), e.getAttribute( "species" ) ) for e in x.getElementsByTagName( "celldesigner:baseProduct" ) ],
                # additional products (same conceptual status as base products, but handled differently in CellDesigner)
                "diss" : [ ( e.getAttribute( "alias" ), e.getAttribute( "product" ) ) for e in x.getElementsByTagName( "celldesigner:productLink" ) ],
            }
            # modifications
            mods = [ ( e.getAttribute( "aliases" ), e.getAttribute( "modifiers" ), e.getAttribute( "type" ) ) for e in x.getElementsByTagName( "celldesigner:modification" ) ]

            if "separate_complexes" in transform:
                # replace complexes with their constituents
                for k, nodes in r.items( ):
                    cnodes = [ ]
                    for n in nodes:
                        if n[ 0 ][ 0 ] == "c":  # select complexes
                            cnodes.extend( [ ( s, s ) for s in sp[ n[ 1 ] ][ "members" ].split( "|" ) ] )
                        else:
                            cnodes.append( n )
                    r[ k ] = cnodes
                cnodes = [ ]
                for i in xrange( len( mods ) ):
                    aliases = mods[ i ][ 0 ].split( "," )
                    species = mods[ i ][ 1 ].split( "," )
                    t = mods[ i ][ 2 ]
                    for n in zip( aliases, species, [ t ] * len( aliases ) ):
                        if n[ 0 ][ 0 ] == "c":  # select complexes
                            cnodes += [ ( n[ 0 ], s, t ) for s in sp[ n[ 1 ] ][ "members" ].split( "|" ) ]
                        else:
                            cnodes.append( n )
                mods = cnodes

            # keep original participants available in edge attributes
            d.update(
                reactants = "|" . join( [ s for s in set( zip( *( r[ "srcs" ] + r[ "adds" ] ) )[ 1 ] ) ] ),
                products = "|" . join( [ t for t in set( zip( *( r[ "tgts" ] + r[ "diss" ] ) )[ 1 ] ) ] ),
                modifiers = "|" . join( [ m[ 1 ] for m in mods ] ),
            )

            # add edges
            for t in r[ "tgts" ] + r[ "diss" ]:
                if t[ 0 ] not in ref2:
                    continue
                # link reagents
                for s in r[ "srcs" ] + r[ "adds" ]:
                    if s[ 0 ] not in ref2:
                        continue
                    try:
                        present = self.graph.edge[ ref2[ s[ 0 ] ] ][ ref2[ t[ 0 ] ] ]
                    except KeyError:
                        present = dict( )
                    type_present = False
                    for e, att in present.iteritems( ):
                        if d[ "type" ] == att[ "type" ]:
                            type_present = True
                            try:
                                revised_reactions[ att[ "r_id" ] ].add( d[ "r_id" ] )
                            except KeyError:
                                revised_reactions[ att[ "r_id" ] ] = { att[ "r_id" ], d[ "r_id" ] }
                    if not type_present:
                        d[ "instance" ] = "{} -D {}" . format( s[ 0 ], t[ 0 ] )
                        self.graph.add_edge( ref2[ s[ 0 ] ], ref2[ t[ 0 ] ], d[ "instance" ], **d )

                # link modifiers
                for m in mods:
                    md = deepcopy( d )  # necessary to copy nested dicts by value instead of by references (as .copy( ) does)
                    md[ "type" ] = m[ 2 ]
                    if md[ "type" ][ -10: ] == "INHIBITION":  # cover INHIBITION and UNKNOWN_INHIBITION
                        md[ "graphics" ].update( target_arrow = 15 )
                        md[ "polarity" ] = -1
                        react = "{} -I {}"
                        if "recolor_edges" in transform:
                            md[ "graphics" ].update( fill = "#ff0000" )
                    else:
                        md[ "graphics" ].update( target_arrow = 12 )
                        md[ "polarity" ] = 1
                        react = "{} -o {}"
                        if "recolor_edges" in transform:
                            md[ "graphics" ].update( fill = "#0000ff" )
                    for a in m[ 0 ].split( "," ):
                        if a not in ref2:
                            continue
                        try:
                            present = self.graph.edge[ ref2[ a ] ][ ref2[ t[ 0 ] ] ]
                        except KeyError:
                            present = dict( )
                        type_present = False
                        for e, att in present.iteritems( ):
                            if md[ "type" ] == att[ "type" ]:
                                type_present = True
                                try:
                                    revised_reactions[ att[ "r_id" ] ].add( d[ "r_id" ] )
                                except KeyError:
                                    revised_reactions[ att[ "r_id" ] ] = { att[ "r_id" ], d[ "r_id" ] }
                        if not type_present:
                            md[ "instance" ] = react.format( a, t[ 0 ] )
                            self.graph.add_edge( ref2[ a ], ref2[ t[ 0 ] ], md[ "instance" ], **md )

        # revise reactions as collected in the previous step
        for r, rset in revised_reactions.iteritems( ):
            revised_reactions[ r ] = "+" . join( sorted( rset, key = lambda x: int( x[ 2: ] ) ) )
        for n1, n2, e in self.graph.edges_iter( data = True ):
            try:
                e[ "r_id" ] = revised_reactions[ e[ "r_id" ] ]
            except KeyError:
                pass

        # produce pseudo reactions if subunits are to be linked
        if "link_subunits" in transform:
            sublink_d = {
                    "r_id" : "psre{}",
                    "metaid" : "",
                    "type" : "COMPLEX_ASSOCIATION",  # newly made up to distinguish it from HETERODIMER_ASSOCIATION
                    "annotation" : "",
                    "graphics" : {
                        "type" : "line",
                        "width" : 4,
                        "fill" : "#880088",
                        "target_arrow" : 6,
                    },
                    "instance" : "",
                    "modifiers" : "",
                    "polarity" : 0,
            }
            # pick out complexes and find neutral edges targeting them - these are the existing association reactions
            nodes = {  k for k, v in self.graph.nodes_iter( data = True ) if v[ "class" ] == "COMPLEX"  }
            linked = dict( [ ( k, set( ) ) for k in nodes ] )
            for n1, n2, e in self.graph.edges_iter( data = True ):
                if n2 in nodes and e[ "polarity" ] == 0:
                    linked[ n2 ].add( self.graph.node[ n1 ][ "reference" ] )
            # check each complex for nodes that need to be added
            i = 1
            increment = 0
            for n in nodes:
                csp = self.graph.node[ n ]
                i += increment  # increment is only 1 if the number i was used in the previous iteration
                increment = 0  # reset each iteration
                sublink_d.update( r_id = "psre" + str( i ), products = csp[ "tag" ] )
                for co in csp[ "members" ].split( "|" ):
                    # do not add edge if the species reference (i.e. basic protein, for example) is already connected
                    if self.graph.node[ ref2[ co ] ][ "reference" ] in linked[ n ]:
                        continue
                    # otherwise add new edge from subunit to complex
                    increment = 1
                    sublink_d.update( instance = "{} -- {}" . format( ref2[ co ], csp[ "alias" ] ), reactants = csp[ "members" ] )
                    self.graph.add_edge( ref2[ co ], ref2[ csp[ "alias" ] ], sublink_d[ "instance" ], **sublink_d.copy( ) )
            self.graph.graph[ "psre" ] = i
        else:
            self.graph.graph[ "psre" ] = 0

        # bridge nodes to be removed from the network
        if "bridge_nodes" in transform and transform[ "bridge_nodes" ]:
            for n in [ n for n in self.graph.nodes_iter( data = True ) if n[ 1 ][ "class" ] in transform[ "bridge_nodes" ] ]:
                self.bridge_node( n )

        # remove specified reactions, or bridge the predecessor node if this removes if this removes its last outgoing edge
        if "bridge_reactions" in transform and transform[ "bridge_reactions" ]:
            for e in self.graph.edges_iter( data = True ):  # iterate over all edges
                if e[ 2 ][ "type" ] in transform[ "bridge_reactions" ]:  # if edges of this type are to be removed
                    self.graph.remove_edge( e[ 0 ], e[ 1 ], e[ 2 ][ "instance" ] )  # NOTE: the specfic key of this edge is not available even with data = True

        # replace "superfamily" nodes with their specified members
        if "replace_nodes" in transform:
            # get nodes environment, i.e. incoming and outgoing edges
            incoming = dict( [ ( k, list( ) ) for k in transform[ "replace_nodes" ] ] )
            outgoing = deepcopy( incoming )
            for e in self.graph.edges_iter( data = True ):
                if e[ 1 ] in incoming:
                    incoming[ e[ 1 ] ].append( e )
                elif e[ 0 ] in outgoing:  # elif prevents self-interactions from being added twice
                    outgoing[ e[ 0 ] ].append( e )
            # add replacement nodes with the same environment
            for node, replacements in transform[ "replace_nodes" ].items( ):
                if node not in self.graph:
                    continue
                for r in replacements:
                    try:
                        rnode = self.graph.node[ r ]
                    except KeyError:
                        d = self.graph.node[ node ].copy( )
                        d.update( tag = r, label = r, name = r[ :-4 ], metaid = node )
                        self.graph.add_node( r, **d )
                        rnode = self.graph.node[ r ]
                    for e in incoming[ node ]:
                        # check if this was a self-loop
                        if e[ 0 ] == node:
                            self.graph.add_edge( r, r, e[ 2 ][ "instance" ], **e[ 2 ] )
                        else:
                            self.graph.add_edge( e[ 0 ], r, e[ 2 ][ "instance" ], **e[ 2 ] )
                    for e in outgoing[ node ]:
                        # no need to check for self-loop because they can only appear in incoming (see above)
                        self.graph.add_edge( r, e[ 1 ], e[ 2 ][ "instance" ], **e[ 2 ] )

            # remove nodes - only after all other operations to prevent nodes from going missing in the middle
            for node in transform[ "replace_nodes" ]:
                if node in self.graph:
                    self.graph.remove_node( node )


    ### ANNOTATION
    def annotate_graph( self ):
        """update nodes with information from the database"""

        if not self.graph:
            self._spill( "Cannot annotate an empty graph." )
            return

        # check each potential annotation for occurrence in the attributes
        refs = dict( )  # to store the mappings per identifier type
        xrefs = dict( self._db._sql( "SELECT `namespace`, `x_id` FROM `Actor_xrefs` UNION SELECT `URN`, `x_id` FROM `Actor_xrefs`" ) )
        xrefs.update( [ ( self._gml_invalid.sub( "_", k ), v ) for k, v in xrefs.iteritems( ) ] )

        # collect which node attributes are present
        attribs = set( )
        for n in self.graph.nodes_iter( data = True ):
            attribs |= set( n[ 1 ] )
        attribs &= set( xrefs )  # reduce by intersecting with available cross-references
        attribs -= { "synonym" }  # remove highly variable identifiers
        if not attribs:
            self._spill( "No node attribues that match any available MIRIAM namespace were found, skipping annotation." )
            return

        # iterate over present node attributes that match one of the MIRIAM namespaces
        for t in attribs:
            query_input = {  str( n[ 1 ][ t ] ) for n in self.graph.nodes_iter( data = True ) if t in n[ 1 ] } - { "" }
            # make sure to convert u'' to normal strings above because there seems to be no sql escape conversion for utf-8 strings
            ref_query = "SELECT `Alias`, `ref` FROM `Actor_aliases` WHERE ref < 1e9 AND type = {} AND Alias IN ({})" . format( xrefs[ t ], self._db._sqllist( query_input ) )
            ref_result = self._db._sql( ref_query )
            if type( ref_result ) == int:
                self._spill( "No recognizable annotation found for identifiers in attribute {}." . format( t ) )
                continue
            refs[ t ] = dict( ref_result )

        # assemble mapped identifiers and requery with them to get mappings
        all_refs = {  v for k in refs for v in refs[ k ].values( )  }
        if not all_refs:
            self._spill( "No recognized annotation found in node attributes, not going to annotate." )
            return

        aliases = defaultdict( lambda: defaultdict( set ) )
        for t in [ self._gml_invalid.sub( "_", x ) for x in ( "uniprot", "ncbigene", "ensembl", "hgnc.symbol", "mirbase.mature" ) ]:
            query = "SELECT Alias, ref FROM Actor_aliases WHERE type = {} AND ref IN ({})" . format( xrefs[ t ], self._db._sqllist( all_refs ) )
            if t == "ensembl":
                query += " AND SUBSTR( Alias, 1, 3 ) = \"ENS\""  # skip old ENSEMBL identifiers
            res = self._db._sql( query )
            for r in res:
                aliases[ r[ 1 ] ][ t ].add( r[ 0 ].strip( ) )

        # update each node with potentially novel annotation
        ref_node = dict( )  # collects references for the subsequent step (complex annotation)
        for n in self.graph.nodes_iter( ):
            node = self.graph.node[ n ]
            try:
                ref_node[ node[ "reference" ] ] = n
            except KeyError:
                pass
            ids = set( )
            for t in refs:
                try:
                    ids.add( refs[ t ][ node[ t ] ] )
                except KeyError:
                    pass

            if len( ids ) > 1:
                self._spill( "Ambiguity detected in mapping for node {}, skipping annotation." . format( n ) )
            elif len( ids ) == 0:
                self._spill( "No regular identifier found for node {}, skipping annotation." . format( n ) )
            else:
                d = dict( )
                e = list( ids )[ 0 ]
                for k, v in aliases[ e ].iteritems( ):
                    if not v:
                        continue
                    v = sorted( v )
                    if len( v ) > 1:
                        d[ k + "_all" ] = "|" . join( v )
                    else:
                        d[ k + "_all" ] = ""
                    if k not in node or not node[ k ]:  # do not overwrite existing annotation
                        d[ k ] = v[ 0 ]
                node.update( d )

                # add this node to the internal symbol mapper dictionary
                for ann in ( "hgnc_symbol", "mirbase_mature" ):
                    if aliases[ e ][ ann ]:
                        symbol = list( aliases[ e ][ ann ] )[ 0 ].lower( )
                        if "_pho" in n:
                            self._symbol_mapper[ "node" ][ symbol ][ "phospho" ].append( n )
                        if "|active" in n:
                            self._symbol_mapper[ "node" ][ symbol ][ "active" ].append( n )
                        else:
                            self._symbol_mapper[ "node" ][ symbol ][ "base" ].append( n )

        # provide complex annotation from annotated nodes
        for n in self.graph.nodes_iter( ):
            node = self.graph.node[ n ]
            if node[ "class" ] != "COMPLEX":  # skip non-complexes
                continue
            elif "subunits" in node:  # skip complexes that have already been processed (as evidenced by the presence of the "subunits" attribute)
                continue
            subunits = set( )
            for m in node[ "members" ].split( "|" ):
                try:
                    subunits.add( self.graph.node[ self._representative[ "spe" ][ m ] ][ "hgnc_symbol" ] )
                except KeyError:
                    continue
            subunits = sorted( subunits )
            node.update( subunits = ":" . join( subunits ) )#, members = "|" . join( subunits ) )

        self._resort_mapper( )

        # update edges with Interactiontype (for Expressence)
        for e in self.graph.edges_iter( data = True ):
            if e[ 2 ][ "polarity" ] == 0:
                e[ 2 ][ "Interactiontype" ] = "interaction"
            elif e[ 2 ][ "polarity" ] == 1:
                e[ 2 ][ "Interactiontype" ] = "stimulation"
            else:
                e[ 2 ][ "Interactiontype" ] = "inhibition"


    ### ANNOTATION FROM EXTERNAL FILE
    def add_annotation( self, files, sep = "\t" ):
        """parse annotation from CSV-like files"""

        if type( files ) not in ( tuple, list, set, dict ):
            self._alert( "Not an iterable: {}" . format( repr( files ) ) )
            return

        for f in files:
            try:
                stream = open( f )
            except IOError:
                self._spill( "Unable to read from file {}." . format( f ) )
                continue

            # header processing
            header = [ self._gml_invalid.sub( "_", h ) for h in stream.readline( ).rstrip( ).split( sep ) ]  # remove newline before splitting
            # prepare the mapping
            if header[ 0 ]:  # the header of the first column names the attribute to map to
                # this is a (potential) one-to-many mapping, so we use a set in the dictionary
                keys = dict( )
                for n in self.graph.nodes_iter( data = True ):
                    try:
                        attrib = n[ 1 ][ header[ 0 ] ]
                    except KeyError:
                        continue
                    else:
                        try:
                            keys[ attrib ].add( n[ 0 ] )
                        except KeyError:
                            keys[ attrib ] = { n[ 0 ] }
            else:  # if the first column's header is empty, we use the internal node id
                # this is always a one-to-one mapping, but it needs to obey the same dictionary structure as the attribute mapping above
                keys = dict( [ ( n, { n } ) for n in self.graph.nodes_iter( ) ] )

            # read file and update nodes
            for line in stream:
                line = line.split( sep )
                try:
                    nodes = keys[ line[ 0 ] ]
                except KeyError:
                    self._spill( "Non-existent node in annotation request: {}." . format( line[ 0 ] ) )
                else:
                    d = dict( [ ( header[ i ], line[ i ].strip( ) ) for i in xrange( 1, len( header ) ) if line[ i ].strip( ) ] )
                    for n in nodes:
                        self.graph.node[ n ].update( d )

            stream.close( )


    ### EXPRESSION VALUES
    def add_exp( self, files, tags = None, impute_complexes = "median-razor", indicate = .05 ):
        """parse expression values from the mentioned files"""

        if type( files ) not in ( tuple, list, set, dict ):
            self._alert( "Not an iterable: {}" . format( repr( files ) ) )
            return

        for i in xrange( len( files ) ):
            f = files[ i ]
            try:
                stream = open( f )
            except IOError:
                self._spill( "Unable to read from file {}." . format( f ) )
                continue

            try:
                tag = tags[ i ]
            except ( TypeError, IndexError ):
                tag = os.path.basename( f ).split( "." )[ 0 ]

            visited = set( )  # keeps track of nodes that have received values from this file (to prevent overwriting by less significant probes)

            # header processing
            stream.readline( )  # discard comment on dataset origin
            header = stream.readline( ).strip( ).split( "\t" )  # get first line
            exp_start = 7  # first column (counting from 0) of expression values
            l = ( len( header ) - exp_start ) / 2  # everything past the first seven columns (gene symbol + limma output) is expression-related
            # check if we can get the sample number from the mean headers (expected format is "HEADER|3")
            for j in xrange( exp_start, exp_start + l ):
                t = header[ j ].split( "|" )
                if len( t ) == 2:
                    try:
                        self.graph.graph[ "_" . join( [ tag, t[ 0 ].split( "_" )[ 0 ], "samples" ] ) ] = int( t[ 1 ] )
                    except ValueError:
                        pass
                header[ j ] = t[ 0 ]
            # clean headers
            header = [ self._gml_invalid.sub( "_", "{}_{}" . format( tag, t ) ) for t in header ]
            logFC = [ h for h in header if "logFC" in h ][ 0 ]

            # read data
            for line in stream:
                line = line.strip( ).split( "\t" )
                symbol = line[ 0 ].lower( )
                if symbol not in self._symbol_mapper[ "node" ]:
                    continue
                elif symbol not in visited:
                    visited.add( symbol )
                    nodes = list( set.union( *[ set( n ) for n in self._symbol_mapper[ "node" ][ symbol ].itervalues( ) ] ) )
                    # add expression values to the fetched nodes
                    exprs = [ ( header[ j ], float( line[ j ] ) ) for j in [ 1, 5 ] + range( 7, len( line ) ) ]
                    self._symbol_mapper[ "expr" ][ symbol ].update( exprs )
                    for n in nodes:
                        self.graph.node[ n ].update( exprs )
                        #self.graph.node[ n ][ logFC.replace( "_logFC", "_FC" ) ] = 2 ** self.graph.node[ n ][ logFC ]

            ### impute values for complexes
            # general preparation
            if impute_complexes:
                # subselect relevant headers
                means = [ h.split( "|" )[ 0 ] for h in header[ 7:(7+l) ] ]
                variances = header[ (7+l): ]

                # specific preparation
                if impute_complexes == "median-razor":  # complex expression = median + differences of subunit from median
                    mod = dict( )  # will hold median values of the sample classes
                    for h in means:
                        #mod[ h ] = median( [ self.graph.node[ n ][ h ] for n in self.graph.nodes_iter( ) if h in self.graph.node[ n ] ] )
                        mod[ h ] = median( [ v[ h ] for v in self._symbol_mapper[ "expr" ].itervalues( ) ] )
                    cimp = lambda iterable, median:  max( 0, sum( iterable ) - ( len( iterable ) - 1 ) * median )
                elif impute_complexes[ :7 ] == "minimum":  # complex expression = minimum of valid subunit expressions
                    mod = dict( [ ( h, 0 ) for h in means ] )  # dummy variable
                    cimp = lambda iterable, waste:  min( set( iterable ) - { 0 } )

                # execution
                for n, att in self.graph.nodes_iter( data = True ):
                    if att[ "class" ] == "COMPLEX":
                        members = [ s.lower( ) for s in att[ "subunits" ].split( ":" ) ]
                        memdicts = [ self._symbol_mapper[ "expr" ][ m ] for m in members ]
                        node = self.graph.node[ n ]
                        fail = False
                        for k in xrange( len( means ) ):
                            h = means[ k ]
                            expr = [ d[ h ] for d in memdicts ]
                            try:
                                i = expr.index( cimp( expr, mod[ h ] ) )
                            except ValueError:  # catches empty set for "minimum"
                                fail = True
                            else:
                                node[ h ] = expr[ i ]
                                node[ variances[ k ] ] = memdicts[ i ][ variances[ k ] ]
                        if not fail:
                            node[ tag + "_logFC" ] = node[ means[ 0 ] ] - node[ means[ 1 ] ]  # NOTE: This will only work for two-class comparisons!
                            #node[ tag + "_FC" ] = 2 ** node[ tag + "_logFC" ]

            if type( indicate ) in ( float, int ):
                if indicate > 0:  # use adjusted p-value
                    try:
                        h = [ h for h in header if "adj_p_val" in h.lower( ) ][ 0 ]
                    except IndexError:
                        continue
                    op = oplt
                else:  # use logFC
                    h = logFC
                    op = opgt
                ind = abs( indicate )
                for n in self.graph.nodes_iter( ):
                    if h in self.graph.node[ n ]:
                        self.graph.node[ n ][ tag + "_indicator" ] = int( op( abs( self.graph.node[ n ][ h ] ), ind ) )


    ### ADDITIONAL INTERACTIONS
    def add_interactions( self, files, criteria = dict( srcs = [ "match:#pro$", "match:[|]active" ], tgts = [ "match:#pro$", "match:_pho[|#]" ] ) ):
        """add interactions from the specified tab-separated CSV files"""

        if type( files ) not in ( tuple, list, set, dict ):
            self._alert( "Not an iterable: {}" . format( repr( files ) ) )
            return

        # prepare criteria for later
        cr = dict( )
        for k, v in criteria.items( ):
            cr[ k ] = list( )
            for e in v:
                p = e.find( ":" )
                rule, value = e[ :p ], e[ p+1: ]
                if rule == "match":
                    cr[ k ].append( re.compile( value ) )

        for f in files:
            try:
                stream = open( f )
            except IOError:
                self._spill( "Unable to read from file {}, skipping." . format( f ) )
                continue

            # header processing
            header = [ self._gml_invalid.sub( "_", h ) for h in stream.readline( ).split( ) ]

            # prepare common node attributes
            n_attrib = dict( graphics = dict( outline_width = 1.0, fill = "#FF00FF", x = 0, y = 0 ) )

            for line in stream:
                line = line.split( "\t" )
                # try to recognize numbers (to provide Cytoscape with the correct type)
                for i in xrange( len( line ) ):
                    try:
                        line[ i ] = int( line[ i ] )
                    except ValueError:
                        try:
                            line[ i ] = float( line[ i ] )
                        except ValueError:
                            pass

                # search for existing nodes to stand in for the new interactors
                nsrcs = [ set( ) ]  # in element n, store set of nodes that have matched all criteria up to n
                ntgts = [ set( ) ]  # dito
                for n in self.graph.nodes_iter( data = True ):
                    if header[ 0 ] in n[ 1 ] and n[ 1 ][ header[ 0 ] ].lower( ) == line[ 0 ].lower( ):
                        nsrcs[ 0 ].add( n[ 0 ] )
                    if header[ 2 ] in n[ 1 ] and n[ 1 ][ header[ 2 ] ].lower( ) == line[ 2 ].lower( ):
                        ntgts[ 0 ].add( n[ 0 ] )

                # select (and possibly add) a source node
                s = None
                if len( nsrcs[ 0 ] ) == 1:
                    s = nsrcs[ 0 ].pop( )
                if len( nsrcs[ 0 ] ) > 1:  # if we have candidates, we single them out according to our criteria
                    i = 0
                    for test in cr[ "srcs" ]:
                        i += 1
                        nsrcs.append( set( ) )
                        for src in nsrcs[ i - 1 ]:
                            if test.search( src ):
                                nsrcs[ i ].add( src )
                        if not nsrcs[ i ]:
                            break
                    while not s and i > -1:
                        try:
                            s = sorted( nsrcs[ i ] )[ 0 ]
                        except IndexError:
                            i -= 1
                # in case of no success above, we need to add the node
                if s == None:
                    s = line[ 0 ].capitalize( )
                    tag = s + "#pro"
                    s_attrib = deepcopy( n_attrib )
                    s_attrib[ header[ 0 ] ] = line[ 0 ]
                    s_attrib.update( name = s, label = tag, tag = tag, homodimer = 1 )
                    s_attrib[ "class" ] = "PROTEIN"
                    self.graph.add_node( tag, **s_attrib )
                    if header[ 0 ] in ( "hgnc_symbol", "mirbase_mature" ):
                        self._symbol_mapper[ "node" ][ s.lower( ) ][ "active" ].append( tag )
                        self._symbol_mapper[ "node" ][ s.lower( ) ][ "base" ].append( tag )
                    s = tag

                # select (and possibly add) a source node
                t = None
                if len( ntgts[ 0 ] ) == 1:
                    t = ntgts[ 0 ].pop( )
                elif len( ntgts[ 0 ] ) > 1:  # if we have candidates, we single them out according to our criteria
                    i = -1
                    for test in cr[ "tgts" ]:
                        i += 1
                        ntgts.append( set( ) )
                        for tgt in ntgts[ i - 1 ]:
                            if test.search( tgt ):
                                ntgts[ i ].add( tgt )
                        if not ntgts[ i ]:
                            break
                    while not t and i > -1:
                        try:
                            t = sorted( ntgts[ i ] )[ 0 ]
                        except IndexError:
                            i -= 1
                # in case of no success above, we need to add the node
                if t == None:
                    t = line[ 2 ].capitalize( )
                    tag = t + "#pro"
                    t_attrib = deepcopy( n_attrib )
                    t_attrib[ header[ 2 ] ] = line[ 2 ]
                    t_attrib.update( name = t, label = tag, tag = tag, homodimer = 1 )
                    t_attrib[ "class" ] = "PROTEIN"
                    self.graph.add_node( tag, **t_attrib )
                    if header[ 2 ] in ( "hgnc_symbol", "mirbase_mature" ):
                        self._symbol_mapper[ "node" ][ t.lower( ) ][ "base" ].append( tag )
                        self._symbol_mapper[ "node" ][ t.lower( ) ][ "phospho" ].append( tag )
                    t = tag

                # select edge attributes based on interaction
                if line[ 1 ] > 0:
                    e_attr = dict( Interactiontype = "stimulation", graphics = dict( fill = "#0000FF", target_arrow = 12, width = 4, type = "line" ) )
                elif line[ 1 ] < 0:
                    e_attr = dict( Interactiontype = "interaction", graphics = dict( fill = "#000000", target_arrow = 6, width = 4, type = "line" ) )
                else:
                    e_attr = dict( Interactiontype = "inhibition", graphics = dict( fill = "#FF0000", target_arrow = 15, width = 4, type = "line" ) )

                e_attr.update( [ ( header[ j ], line[ j ] ) for j in ( 1, 3, 4 ) ] )
                self.graph.graph[ "psre" ] += 1
                e_attr.update( r_id = "psre{}" . format( self.graph.graph[ "psre" ] ), type = "CURATED", instance = "{} -o {}" . format( s, t ) )
                self.graph.add_edge( s, t, e_attr[ "instance" ], **e_attr )

            stream.close( )


############################################################
#### PRIVATE                                            ####
############################################################

    def _reset( self ):
        """reset internal storages"""

        NetworkCreator._reset( self )
        self.xml_object = None
        self._representative = dict( spe = dict( ), ref = dict( ), ali = dict( ) )


    ### HELPER FUNCTION
    def bridge_node( self, node ):
        """connect predecessors and successors of the node, then remove it"""

        neighbors = {  "incoming" : set( ), "outgoing" : set( )  }
        # get successors directly from edge relations
        if node[ 0 ] in self.graph.edge:
            neighbors[ "outgoing" ] = set( self.graph.edge[ node[ 0 ] ] ) & set( self.graph ) - { node[ 0 ] }  # check for removed nodes and self-loops
        # iterator over edges to identify predecessors
        for e in self.graph.edges_iter( ):
            if e[ 1 ] == node[ 0 ] and e[ 1 ] in self.graph and e[ 0 ] != node[ 0 ]:  # check for removed nodes and self-loops
                neighbors[ "incoming" ].add( e[ 0 ] )

        ga = {  "source_arrow" : 0, "target_arrow" : 12, "fill" : "#ff00ff", "width" : 4, "type" : "line"  }  # activating edge
        gn = {  "source_arrow" : 0, "target_arrow" : 6, "fill" : "#ff00ff", "width" : 4, "type" : "line"  }  # neutral edge
        gi = {  "source_arrow" : 0, "target_arrow" : 15, "fill" : "#ff00ff", "width" : 4, "type" : "line"  }  # inhibitory edge
        # iterate over the neighbors and connect them among each other with the appropriate edges.
        for pred in neighbors[ "incoming" ]:
            pred_edges = self.graph.edge[ pred ][ node[ 0 ] ].copy( )  # dict of edges between predecessor and node

            for succ in neighbors[ "outgoing" ] - neighbors[ "incoming" ]:  # in reversible reactions (PIP2<>PIP3 and others), avoid briding the inverse leg
                succ_edges = self.graph.edge[ node[ 0 ] ][ succ ].copy( )  # dict of edges between node and successor

                for p, pv in pred_edges.items( ):
                    for s, sv in succ_edges.items( ):
                        d = dict( [ ( k, "{}+{}" . format( pv[ k ], sv[ k ] ) ) for k in set( pv ) - set( [ "graphics" ] ) ] )
                        d[ "type" ] = "BRIDGED:" + d[ "type" ].replace( "BRIDGED:", "" )  # remove existing BRIDGED:s before prefixing one
                        # calculate how many "inhibitory" edges are between the two nodes: 1 makes the new edge inhibitory, 0 or 2 make it activating
                        if pv[ "polarity" ] == sv[ "polarity" ]:
                            if pv[ "polarity" ] == 0:
                                d.update( polarity = 0, graphics = gn )
                            else:
                                d.update( polarity = 1, graphics = ga )
                        elif -1 in ( pv[ "polarity" ], sv[ "polarity" ] ):
                            d.update( polarity = -1, graphics = gi )
                        else:
                            d.update( polarity = 1, graphics = ga )
                        self.graph.add_edge( pred, succ, "{}+{}" . format( p, s ), **d )

        # remove the matching node - automatically removes all edges connected to that node
        self.graph.remove_node( node[ 0 ] )



