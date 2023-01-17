#!/usr/bin/env bash
# Copyright [2021] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

## A script to help download and process data from OrthoDB [See: https://www.orthodb.org/orthodb_userguide.html#contact]
## Author: Lahcen Campbell [lcampbell@ebi.ac.uk]
## version 1.4

## Update 1.4: -Improved handling of main download URL and orthoDB version needed for OrthoDB master file/data retreval
##             -Implement embedding of OrthoDB version into protein file names

# Main input variable
TAXID_CLADE=$1
CWD=`readlink -f $PWD`
PERL_SCRIPTS_DIR="../support_scripts_perl"

## IMPORTANT URL - Which could be changed in a future OrthoDB update.... ##
ORTHODB_FILE_URL="https://data.orthodb.org/download"
# [ReadMe = $ORTHODB_FILE_URL/README.txt]

if [ -z $TAXID_CLADE ]; then
echo -e -n 'Taxon ID OR Clade name required. Exiting...\n\nUsage: sh download_orthodb_proteinset.sh <TaxonID -OR- Clade Name>\n'
echo -e -n 'E.g:\nsh Download_OrthoDB_ProtSet.sh mollusca\nOR\n'
echo 'sh Download_OrthoDB_ProtSet.sh 6447'
exit 1
fi

# Function for testing taxonID vs clade name input
is_int () { test "$@" -eq "$@" 2> /dev/null; }

## Get the lastest information related to the current version set of *.tab.gz files hosted on OrthoDB:
wget -q $ORTHODB_FILE_URL -O OrthoDB_Download.html
ODB_LEVEL2SPECIES=`grep -e '_level2species.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+_level2species.tab.gz/){print $&."\n";};' | head -n 1`
ODB_LEVELS=`grep -e '_levels.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+_levels.tab.gz/){print $&."\n";};' | head -n 1`
ODB_VERSION=`grep -e '_all_fasta.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+/ ){print $&."\n";};' | head -n 1`
echo "Using OrthoDB Version: $ODB_VERSION"

# Test for presence of non fasta files from OrthoDB. Used to gain clade/species information.
if [[ ! -f ${CWD}/$ODB_LEVEL2SPECIES ]] && [[ ! -f ${CWD}/$ODB_LEVELS ]]; then
    echo "## Downloading OrthoDB taxonomy to Ortho master files...."
    for ORTHFILE in $ODB_LEVEL2SPECIES $ODB_LEVELS
    do
        wget -q ${ORTHODB_FILE_URL}/$ORTHFILE
    done
fi

## Main processing begings here onwards...
if [[ -v "$TAXON_CLADE" ]]; then
    echo "TaxonID or Clade name not provided. EXITING"
    exit 0
elif is_int "$TAXID_CLADE"; then
    echo "## Processing on Taxon ID: '$TAXID_CLADE'"
    TAXON_ID=$TAXID_CLADE
    CLADE_NAME=`zcat $ODB_LEVELS | grep -w -i -e "$TAXID_CLADE" | cut -f2`
    CLADE_TID_LIST=`zcat $ODB_LEVEL2SPECIES | grep -w -e "$TAXID_CLADE" | cut -f4 | awk 'BEGIN { FS="," } { print $NF }' | sed 's/}//' | tr "\n" "," | sed 's/,$//g'`
    if [[ -z $CLADE_TID_LIST ]]; then 
        echo "!!! Taxon ID: $TAXID_CLADE not defined within OrthoDB. See $ODB_LEVELS"
        exit 1
    else
        #Obtain taxon information from uniprot
        echo $CLADE_TID_LIST | tr "," "\n" | xargs -n 1 -I XXX wget -q 'https://rest.uniprot.org/taxonomy/XXX.tsv' -O ->> ${CLADE_NAME}.comb.uniprot.tmp
        head -n 1 ${CLADE_NAME}.comb.uniprot.tmp >> ${CLADE_NAME}.orthodb.uniprot.tsv
        grep -E "^[0-9]" ${CLADE_NAME}.comb.uniprot.tmp >> ${CLADE_NAME}.orthodb.uniprot.tsv
        rm ${CLADE_NAME}.comb.uniprot.tmp
    fi
else
    echo "## Processing on CLADE name '$TAXID_CLADE'"
    CLADE_NAME=$TAXID_CLADE
    TAXON_ID=`zcat $ODB_LEVELS | grep -w -i -e "$TAXID_CLADE" | cut -f1`
    CLADE_TID_LIST=`zcat $ODB_LEVEL2SPECIES | grep -w -e "$TAXON_ID" | cut -f4 | awk 'BEGIN { FS="," } { print $NF }' | sed 's/}//' | tr "\n" "," | sed 's/,$//g'`
    if [[ -z $TAXON_ID ]] || [[ -z $CLADE_TID_LIST ]]; then 
        echo "!!! Clade name: '$TAXID_CLADE' not located within OrthoDB. See $ODB_LEVELS"
        exit 1
    else
        #Obtain taxon information from uniprot
        echo $CLADE_TID_LIST | tr "," "\n" | xargs -n 1 -I XXX wget -q 'https://rest.uniprot.org/taxonomy/XXX.tsv' -O ->> ${CLADE_NAME}.comb.uniprot.tmp
        head -n 1 ${CLADE_NAME}.comb.uniprot.tmp >> ${CLADE_NAME}.orthodb.uniprot.tsv
        grep -E "^[0-9]" ${CLADE_NAME}.comb.uniprot.tmp >> ${CLADE_NAME}.orthodb.uniprot.tsv
        rm ${CLADE_NAME}.comb.uniprot.tmp
    fi
fi

# Ensure we dont have a species level taxonID. Taxon ID or clade name can not be a single species.
if [[ -z "$CLADE_NAME" ]]; then
    echo "!!! Clade name using '$TAXID_CLADE' not located. See [ $ODB_LEVEL2SPECIES = 'Correspondence between level ids and organism ids' ]."
    if is_int "$TAXON_ID"; then
    `wget -q https://rest.uniprot.org/taxonomy/${TAXON_ID}.tsv -O Uniprot_taxid${TAXON_ID}.tsv`
    RANK=`tail -n 1 Uniprot_taxid${TAXON_ID}.tsv | cut -f7`
    # Test we are not using species level taxon rank to get clusters
        if [[ $RANK =~ "species" ]]; then
            echo "Taxon ID is at SPECIES level. Don't use species level taxIDs. Try Sub/Infra/Order level ID or higher."
        fi
    fi
    echo "## See File: Uniprot_taxid${TAXON_ID}.tsv"
    exit 1
fi

# Report set of taxon IDs to screen
LOG_CLUSTERS="${CWD}/${CLADE_NAME}_orthodb_download.cluster.log.txt"
TAXON_COUNT=`grep -c -E "^[0-9]" ${CLADE_NAME}.orthodb.uniprot.tsv`
echo -e -n "$CLADE_NAME clade contains set of $TAXON_COUNT Taxon IDs => [ $CLADE_TID_LIST ]\n" | tee -a $LOG_CLUSTERS
echo -e -n "\tSee Taxon information from uniprot in file: ${CLADE_NAME}.orthodb.uniprot.tsv\n\n" | tee -a $LOG_CLUSTERS

## Get the set of OrthoDB clusters based on clade of interest using taxonID info
CLUSTER_IDS_URL="https://v101.orthodb.org//search?query=&level=${TAXON_ID}&species=${TAXON_ID}&universal=0.9&singlecopy=&limit=80000"
# Download cluster json using Taxon ID
wget $CLUSTER_IDS_URL -O ${CLADE_NAME}_orthoDB_clusters.json 2> /dev/null
echo "wget $CLUSTER_IDS_URL -O ${CLADE_NAME}_orthoDB_clusters.json"

## Set output file name for combined clusters text file 
LINEAR_CLUSTERS="${CLADE_NAME}_orthoDB_clusters.linear.txt"
jq '.data[]' ${CLADE_NAME}_orthoDB_clusters.json > ${CWD}/$LINEAR_CLUSTERS
sed -i 's/"//g' ${CWD}/$LINEAR_CLUSTERS # File used to parse clusters and download them individually
CLUSTER_COUNT=`jq '.count' ${CLADE_NAME}_orthoDB_clusters.json`
echo -e -n "\n*** Retrieved a total of $CLUSTER_COUNT clusters for $CLADE_NAME ***\n\n"

sleep 3

## Download individual clusters using the set of cluster IDs in XX and the CLADE_TAXID set
ORIG_CLUSTERS_COMB="${CWD}/Combined_${CLADE_NAME}_OrthoDB.orig.fa"
if [ -f ${CWD}/Combined_${CLADE_NAME}_OrthoDB.orig.fa ]; then
    rm ${CWD}/Combined_${CLADE_NAME}_OrthoDB.orig.fa
fi

# Download individual OrthoDB clusters
while read CLUSTER
do
    SINGLE_CLUSTER="${CWD}/sub_${CLUSTER}.fa"
    echo "wget -qq 'https://v101.orthodb.org/fasta?query=level=${TAXON_ID}&id=${CLUSTER}&species=${CLADE_TID_LIST}' -O $SINGLE_CLUSTER" > ${CWD}/OrthoDB_wget_${CLUSTER}.sh
    echo -e -n "Processing OrthoDB cluster: $CLUSTER\n"
    echo -e -n "Running --> OrthoDB_wget_${CLUSTER}.sh\n" >> $LOG_CLUSTERS
    cat OrthoDB_wget_${CLUSTER}.sh >> $LOG_CLUSTERS;
    sh ${CWD}/OrthoDB_wget_${CLUSTER}.sh 2>&1 | tee -a $LOG_CLUSTERS
    cat $SINGLE_CLUSTER >> $ORIG_CLUSTERS_COMB; rm $SINGLE_CLUSTER
    rm  ${CWD}/OrthoDB_wget_${CLUSTER}.sh
done < $LINEAR_CLUSTERS

## Process the Combined cluster DB fasta to sort out headers. Retaining unique orthoDB seqIDs, but removing...
## redundancy and adding a counter when sequences are not unique to single OrthoDB clusters.
ORIG_CLUST_JUSTFILE="Combined_${CLADE_NAME}_OrthoDB.orig.fa"
BASENAME=(`basename $ORIG_CLUST_JUSTFILE .orig.fa`)
FINAL_ORTHO_FASTA="${CWD}/${BASENAME}_final.out.fa"
DEDUP_OUT_TMP="${CWD}/${BASENAME}_no_dups.tmp"
REHEADER_OUT="${CWD}/${BASENAME}_Reheader.fa.tmp"

## Process original cluster file to remove any duplicated sequences prior to downstream processing
echo -e -n "\n## Removing duplicate sequences from [ $ORIG_CLUST_JUSTFILE ]\n"
perl ${PERL_SCRIPTS_DIR}/remove_dup_seqs.pl $ORIG_CLUST_JUSTFILE > $DEDUP_OUT_TMP

## Reduce the set of headers in the OrthoDB fasta file to the unique sequence Identifier.
echo "## Processing deduplicated seq headers, isolating to unique OrthoDB gene ID..."
perl ${PERL_SCRIPTS_DIR}/reheader_orthodb.pl $DEDUP_OUT_TMP $BASENAME

# Create final outputfile based on reheadered orthoDB file
mv $REHEADER_OUT $FINAL_ORTHO_FASTA

## Generate samtools index file
echo "## Creating samtools index of $FINAL_ORTHO_FASTA"
`samtools faidx $FINAL_ORTHO_FASTA`

echo -e -n "\n*** Processing OrthoDB for [ $CLADE_NAME ] completed !! ***\n\n"
rm ${CWD}/*.tmp
echo "## Fasta sequence count in original and final fasta file"
grep -c -e ">" -I ${BASENAME}*.fa

## Print command to rename files:
echo -e -n "\n## Renaming final output files:\n"
perl ${PERL_SCRIPTS_DIR}/Quick_rename.pl ${BASENAME}_final.out. ${CLADE_NAME}_orth${ODB_VERSION}_proteins. prefix

exit 1
