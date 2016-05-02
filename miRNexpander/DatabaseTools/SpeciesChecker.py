#!/usr/bin/python

# this script checks the species for those missing species in the miRNA name

import re

Species_Dict = {  # one of the list elements will constitute the canonical name that will be returned by get_canonical_name ( can be changed, see below )
    "hsa" : [ 9606, 'hsa', 'human', 'homo sapiens', 'h. sapiens', 'homo' ],
    "mmu" : [ 10090, 'mmu', 'mouse', 'mus musculus', 'm. musculus', 'mus', 'murine' ],
    "rno" : [ 10116, 'rno', 'rat', 'rattus norvegicus', 'r. norvegicus', 'rattus' ],
    "dme" : [ 7227, 'dme', 'fruitfly', 'drosophila melanogaster', 'd. melanogaster', 'drosophila' ],
    "gga" : [ 9031, 'gga', 'chicken', 'gallus gallus', 'g. gallus', 'gallus', 'chick' ],
    "dre" : [ 7955, 'dre', 'zebrafish', 'danio rerio', 'd. rerio', 'danio', 'brachydanio rerio' ],
    "cel" : [ 6239, 'cel', 'roundworm', 'caenorhabditis elegans', 'c. elegans', 'caenorhabditis' ],
    "oar" : [ 9940, 'oar', 'sheep', 'ovis aries', 'o. aries', 'ovis' ],
    "ocu" : [ 9986, 'ocu', 'rabbit', 'oryctolagos caniculus', 'o. caniculus', 'oryctolagos' ],
    "bta" : [ 9913, 'bta', 'cattle', 'bos taurus', 'b. taurus', 'bos', 'bos primigenius', 'bos primigenius taurus' ],
    "osa" : [ 4530, 'osa', 'asian rice', 'oryza sativa', 'o. sativa', 'oryza', 'rice' ],
    "ogl" : [ 4538, 'ogl', 'african rice', 'oryza glaberrima', 'o. glaberrima' ],
    "ath" : [ 3702, 'ath', 'thale cress', 'arabidopsis thaliana', 'a. thaliana', 'arabidopsis', 'mouse-ear cress' ]
# consider adding pig/swine, horse, dog, cat...
}

Viruses = re.compile( 'epstein-barr|papilloma' )  # this definitely needs expansion


def get_species_id( species ):
    """return the GenBank Identifier for some common species"""

    species = species.strip( ).lower( )
    result = 1  # (non-sensical) fail-safe if there is no match in the loop
    for species_key in Species_Dict:
        if species in Species_Dict[ species_key ]:
            result = Species_Dict[ species_key ][ 0 ]  # change assignment if you want to return another list element
            break
    return result


def get_scientific_name( species ):
    """return the scientific species name for some common species"""

    species = species.strip( ).lower( )
    result = species  # fail-safe if there is no match in the loop
    for species_key in Species_Dict:
        if species in Species_Dict[ species_key ]:
            result = Species_Dict[ species_key ][ 3 ]  # change assignment if you want to return another list element
            break
    return result


def get_latin_abbr( species ):
    """return a 3-character abbreviation of the given species"""

    species = species.strip( ).lower( )
    result = "xxx"  # fail-safe if there is no match
    if Viruses.search( species ):
        result = "vir"
    else:
        for species_key in Species_Dict:
            if species in Species_Dict[ species_key ]:
                result = species_key  # change assignment if you want to return a list member, e.g. result = Species_Dict[ species_key ][ 0 ]
                break
    return result


