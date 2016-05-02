#!/bin/bash

# This file corrects spelling and structural errors in the output of some databases.
# Comprehensiveness is not guaranteed.

### HTRIdb_data.txt
# tabulator mess
f=./HTRI/HTRIdb_data.txt
sed -r -i".bak" "s/\t+/\t/g" "$f"


### miRecords_version4.csv
# misspellings in some miRNA names
f=./miRecords/miRecords_version4.csv
sed -r -i".bak" "/Homo sapiens/ s/has-/hsa-/
                " "$f"

### r20_miRBase.dat
# misspellings, ambiguity, or unauthoritative variants in taxonomic identifiers
# NOTE: not all ambiguities might be resolved perfectly
f=./miRBase/r20_miRNA.dat
sed -r -i".bak" "/^DE/ s/Mareks disease/Marek's disease/
                    /^DE/ s/Kaposis? sarcoma-associated/Kaposi's sarcoma-associated/
                    /^DE/ s/Mouse gammaherpesvirus 68/murine gammaherpesvirus 68/
                    /^DE/ s/Mouse cytomegalovirus/murine cytomegalovirus/
                    /^DE/ s/Saccharum ssp./sugarcane <hybrid>/
                    /^DE/ s/Herpesvirus saimiri strain A11/Herpesvirus saimiri (strain 11)/
                    /^DE/ s/Citrus clementine/Citrus clementina/" "$f"

