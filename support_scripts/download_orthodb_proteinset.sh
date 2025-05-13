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

## API Docs: https://www.ezlab.org/orthodb_userguide.html
## Author: Lahcen Campbell [lcampbell@ebi.ac.uk]
## Current version: 1.6

## Updated in 1.6 - Updated Rest Endpoint URLs (as of April 2025)
##                - Increased max cluster retrieval to new max 10k (up from 5K)
##                - Added small check for path access to support_scripts dir
##                - Assorted spelling fixes, variable formatting consistency

# Main input variable
TAXID_CLADE=${1^}
ORTHODB_VERSION="v"$2
CWD_TMP=`readlink -f $PWD`
DATE_RUN=$(date "+%d/%m/%Y" | sed 's/\//-/g')
PERL_SCRIPTS_DIR="$CWD_TMP/ensembl-anno/support_scripts_perl"
PROCESSING_LOG="Trace.${DATE_RUN}.log"

## IMPORTANT URL - Which could be changed in a future OrthoDB update.... ##
ORTHODB_FILE_URL="https://data.orthodb.org/current/download"
# [ReadMe = $ORTHODB_FILE_URL/README.txt]

# Define endpoint URL if user specified version:
if [[ "${ORTHODB_VERSION}" =~ v[0-9]+ ]]; then
    ORTHODB_FILE_URL="https://data.orthodb.org/${ORTHODB_VERSION}/download"
fi

## check URL is active:
HTTP_CODE=`curl -sL -w "%{http_code}\n" "$ORTHODB_FILE_URL" -o /dev/null`
if [[ $HTTP_CODE != "200" ]]; then
    echo "URL - > $ORTHODB_FILE_URL is returning HTTP code $HTTP_CODE"
    echo -e -n "Is this expected ?\nNeed to exit..."
    exit
else
    echo -e -n "$ORTHODB_FILE_URL has non-404 status hurray (code:$HTTP_CODE)\nProceeding...\n" | tee $PROCESSING_LOG
fi

# Check before doing any processing if support scripts dir is available
if [[ ! -d $PERL_SCRIPTS_DIR ]]; then
    echo -e -n "Unable to locate and define perl support scripts directory variable '${PERL_SCRIPTS_DIR}'.\n"
    echo -e -n "Please make sure you have this path available, set 'PERL_SCRIPTS_DIR':line35.\nExiting.\n"
    exit
fi

if [ -z $TAXID_CLADE ]; then
    echo -e -n 'Taxon ID OR Clade name required. Exiting...\n\nUsage: sh download_orthodb_proteinset.sh <TaxonID -OR- Clade Name> <Optional: OrthoDB version [e.g 12]>\n'
    echo -e -n 'E.g:\nsh Download_OrthoDB_ProtSet.sh mollusca 12\nOR\nsh Download_OrthoDB_ProtSet.sh mollusca\nOR\n'
    echo 'sh Download_OrthoDB_ProtSet.sh 6447'
    exit 1
else
    CWD=$CWD_TMP/OrthoDB_${DATE_RUN}_${TAXID_CLADE}
    mkdir -p $CWD
    cd $CWD
    echo -e -n "Processing start: " | tee -a $PROCESSING_LOG
    date | tee -a $PROCESSING_LOG
fi

# Function for testing taxonID vs clade name input
is_int () { test "$@" -eq "$@" 2> /dev/null; }

## Get the latest information related to the current version set of *.tab.gz files hosted on OrthoDB:
wget -q $ORTHODB_FILE_URL -O OrthoDB_Download.html
ODB_LEVEL2SPECIES=`grep -e '_level2species.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+_level2species.tab.gz/){print $&."\n";};' | head -n 1`
ODB_LEVELS=`grep -e '_levels.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+_levels.tab.gz/){print $&."\n";};' | head -n 1`
ODB_VERSION=`grep -e '_all_fasta.tab.gz' OrthoDB_Download.html | perl -pe 'if ( $_ =~ m/odb[0-9]+v[0-9]+/ ){print $&."\n";};' | head -n 1`
echo "Using OrthoDB Version: $ODB_VERSION" | tee -a $PROCESSING_LOG

# Test for presence of non fasta files from OrthoDB. Used to gain clade/species information.
if [[ ! -f ${CWD}/$ODB_LEVEL2SPECIES ]] && [[ ! -f ${CWD}/$ODB_LEVELS ]]; then
    echo "## Downloading OrthoDB taxonomy to Ortho master files...." | tee -a $PROCESSING_LOG
    for ORTHFILE in $ODB_LEVEL2SPECIES $ODB_LEVELS
    do
        wget -q ${ORTHODB_FILE_URL}/$ORTHFILE
    done
fi

## Main processing begings here onwards...
if [[ -v "$TAXON_CLADE" ]]; then
    echo "TaxonID or Clade name not provided. EXITING"  | tee -a $PROCESSING_LOG
    exit 0
elif is_int "$TAXID_CLADE"; then
    echo "## Processing on Taxon ID: '$TAXID_CLADE'"  | tee -a $PROCESSING_LOG
    TAXON_ID=$TAXID_CLADE
    CLADE_NAME=`zcat $ODB_LEVELS | grep -w -i -e "$TAXID_CLADE" | cut -f2`
    CLADE_TID_LIST=`zcat $ODB_LEVEL2SPECIES | grep -w -e "$TAXID_CLADE" | cut -f4 | awk 'BEGIN { FS="," } { print $NF }' | sed 's/}//' | tr "\n" "," | sed 's/,$//g'`
    if [[ -z $CLADE_TID_LIST ]]; then 
        echo "!!! Taxon ID: $TAXID_CLADE not defined within OrthoDB. See $ODB_LEVELS" | tee -a $PROCESSING_LOG
        exit 1
    else
        #Obtain taxon information from uniprot
        echo $CLADE_TID_LIST | tr "," "\n" | xargs -n 1 -I XXX wget -q 'https://rest.uniprot.org/taxonomy/XXX.tsv' -O ->> ${CLADE_NAME}.comb.uniprot.tmp
        head -n 1 ${CLADE_NAME}.comb.uniprot.tmp > ${CLADE_NAME}.orthodb.uniprot.tsv
        grep -E "^[0-9]" ${CLADE_NAME}.comb.uniprot.tmp >> ${CLADE_NAME}.orthodb.uniprot.tsv
        rm ${CLADE_NAME}.comb.uniprot.tmp
    fi
else
    echo "## Processing on CLADE name '$TAXID_CLADE'" | tee -a $PROCESSING_LOG
    CLADE_NAME=$TAXID_CLADE
    TAXON_ID=`zcat $ODB_LEVELS | grep -w -i -e "$TAXID_CLADE" | cut -f1`
    CLADE_TID_LIST=`zcat $ODB_LEVEL2SPECIES | grep -w -e "$TAXON_ID" | cut -f4 | awk 'BEGIN { FS="," } { print $NF }' | sed 's/}//' | tr "\n" "," | sed 's/,$//g'`
    if [[ -z $TAXON_ID ]] || [[ -z $CLADE_TID_LIST ]]; then 
        echo "!!! Clade name: '$TAXID_CLADE' not located within OrthoDB. See $ODB_LEVELS" | tee -a $PROCESSING_LOG
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
    echo "!!! Clade name using '$TAXID_CLADE' not located. See [ $ODB_LEVEL2SPECIES = 'Correspondence between level ids and organism ids' ]." | tee -a $PROCESSING_LOG
    if is_int "$TAXON_ID"; then
    `wget -q https://rest.uniprot.org/taxonomy/${TAXON_ID}.tsv -O Uniprot_taxid${TAXON_ID}.tsv`
    RANK=`tail -n 1 Uniprot_taxid${TAXON_ID}.tsv | cut -f7`
    # Test we are not using species level taxon rank to get clusters
        if [[ $RANK =~ "species" ]]; then
            echo "Taxon ID is at SPECIES level. Don't use species level taxIDs. Try Sub/Infra/Order level ID or higher."
        fi
    fi
    echo "## See File: Uniprot_taxid${TAXON_ID}.tsv"  | tee -a $PROCESSING_LOG
    exit 1
fi

# Report set of taxon IDs to screen
LOG_CLUSTERS="${CWD}/clusters_OrthoDB_download.${CLADE_NAME}.log.txt"
TAXON_COUNT=`grep -c -E "^[0-9]" ${CLADE_NAME}.orthodb.uniprot.tsv`
echo -e -n "$CLADE_NAME clade contains set of $TAXON_COUNT Taxon IDs => [ $CLADE_TID_LIST ]\n" | tee -a $PROCESSING_LOG
echo -e -n "\tSee Taxon information from uniprot in file: ${CLADE_NAME}.orthodb.uniprot.tsv\n\n" | tee -a $PROCESSING_LOG

## Get the set of OrthoDB clusters based on clade of interest using taxonID info
CLUSTER_IDS_URL="https://data.orthodb.org/current/search?universal=0.9&singlecopy=&level=${TAXON_ID}&species=${TAXON_ID}&take=10000"

## Download cluster json using Taxon ID
CLUSTER_JSON="${CLADE_NAME}.orthodb_clusters.json"
echo "curl $CLUSTER_IDS_URL -o $CLUSTER_JSON 2> /dev/null" | tee -a $PROCESSING_LOG
curl $CLUSTER_IDS_URL -L -o $CLUSTER_JSON 2> /dev/null

## Set output file name for combined clusters text file 
LINEAR_CLUSTERS="${CWD}/linear_${CLADE_NAME}_orthoDB_clusters.txt"
jq '.data[]' $CLUSTER_JSON | sed 's/"//g' > $LINEAR_CLUSTERS  # File used to parse clusters and download them individually
CLUSTER_COUNT=`jq '.count' ${CWD}/${CLUSTER_JSON}`
echo -e -n "\n*** Retrieved a total of $CLUSTER_COUNT clusters for $CLADE_NAME ***\n\n" | tee -a $PROCESSING_LOG
sleep 3

## Download individual clusters using the set of cluster IDs in XX and the CLADE_TAXID set
ORIG_CLUSTERS_COMB="${CWD}/Combined_${CLADE_NAME}_OrthoDB.orig.fa"
if [ -f $ORIG_CLUSTERS_COMB ]; then
    mv $ORIG_CLUSTERS_COMB ${CWD}/Combined_${CLADE_NAME}_OrthoDB.orig.OLD.fa
fi

# Download individual OrthoDB clusters and combin them into a single fasta file (contains duplicated sequences)
while read CLUSTER_ID
do
    # Cluster Variables
    SINGLE_CLUSTER_FA="${CWD}/sub_${CLUSTER_ID}.fa"
    FASTA_CLUSTER_URL="https://data.orthodb.org/current/fasta?id=${TAXON_ID}&id=${CLUSTER_ID}&species=${CLADE_TID_LIST}"
    CURL_SCRIPT="${CWD}/OrthoDB_curl_${CLUSTER_ID}.sh"
    echo "curl -s \"$FASTA_CLUSTER_URL\" -o $SINGLE_CLUSTER_FA" > $CURL_SCRIPT
    echo -e -n "Processing OrthoDB cluster: $CLUSTER_ID\n" | tee -a $PROCESSING_LOG
    echo -e -n "Running --> $CURL_SCRIPT\n" >> $LOG_CLUSTERS
    cat $CURL_SCRIPT >> $LOG_CLUSTERS
    sh $CURL_SCRIPT 2>&1 | tee -a $LOG_CLUSTERS
    cat $SINGLE_CLUSTER_FA >> $ORIG_CLUSTERS_COMB
    rm $SINGLE_CLUSTER_FA
    rm $CURL_SCRIPT
done < $LINEAR_CLUSTERS

## Process the Combined cluster DB fasta to sort out headers. Retaining unique orthoDB seqIDs, but removing...
## redundancy and adding a counter when sequences are not unique to single OrthoDB clusters.
ORIG_CLUST_JUSTFILE="Combined_${CLADE_NAME}_OrthoDB.orig.fa"
BASENAME=(`basename $ORIG_CLUST_JUSTFILE .orig.fa`)
FINAL_ORTHO_FASTA="${CWD}/${BASENAME}_final.out.fa"
DEDUP_OUT_TMP="${CWD}/${BASENAME}_no_dups.tmp"
REHEADER_OUT="${CWD}/${BASENAME}_Reheader.fa.tmp"

## Process original cluster file to remove any duplicated sequences prior to downstream processing
echo -e -n "\n## Removing duplicate sequences from [ $ORIG_CLUST_JUSTFILE ]\n" | tee -a $PROCESSING_LOG
perl ${PERL_SCRIPTS_DIR}/remove_dup_seqs.pl $ORIG_CLUST_JUSTFILE > $DEDUP_OUT_TMP

## Reduce the set of headers in the OrthoDB fasta file to the unique sequence Identifier.
echo "## Processing deduplicated seq headers, isolating to unique OrthoDB gene ID..." | tee -a $PROCESSING_LOG
perl ${PERL_SCRIPTS_DIR}/reheader_orthodb.pl $DEDUP_OUT_TMP $BASENAME

# Create final output file based on reheadered orthoDB file
mv $REHEADER_OUT $FINAL_ORTHO_FASTA

## Generate samtools index file
echo "## Creating samtools index of $FINAL_ORTHO_FASTA" | tee -a $PROCESSING_LOG
`samtools faidx $FINAL_ORTHO_FASTA`

echo -e -n "\n*** Processing OrthoDB for [ $CLADE_NAME ] completed !! ***\n\n" | tee -a $PROCESSING_LOG
rm ${CWD}/*.tmp
echo "## Fasta sequence count in original combined fasta [Includes Duplicates seqs]" | tee -a $PROCESSING_LOG
echo -e -n "$ORIG_CLUST_JUSTFILE\tNum seqs. --> " | tee -a $PROCESSING_LOG
grep -c -e ">" $ORIG_CLUST_JUSTFILE | tee -a $PROCESSING_LOG

## Print command to rename files:
echo -e -n "\n## Renaming final output files:\n" | tee -a $PROCESSING_LOG
perl ${PERL_SCRIPTS_DIR}/quick_rename.pl ${BASENAME}_final.out. ${CLADE_NAME}_orth${ODB_VERSION}_proteins.uniq. prefix 2>&1> /dev/null

echo "## Final OrthoDB fasta outfile seq count [Unique sequences]" | tee -a $PROCESSING_LOG
echo -e -n "${CLADE_NAME}_orth${ODB_VERSION}_proteins.uniq.fa\tNum seqs. --> " | tee -a $PROCESSING_LOG
grep -c -e ">" ${CLADE_NAME}_orth${ODB_VERSION}_proteins.uniq.fa | tee -a $PROCESSING_LOG

echo -e -n "Processing finished: " | tee -a $PROCESSING_LOG
date | tee -a $PROCESSING_LOG
exit 1
