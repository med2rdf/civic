#!/bin/bash
DOCKER_DIR="docker_build_c"

cp -p tsv2rdf_civic.py $DOCKER_DIR
cp -p templ_civic.ttl $DOCKER_DIR
cp -p templ_civic.ttl.evi $DOCKER_DIR
cp -p templ_civic.ttl.prefix $DOCKER_DIR
cp -p templ_civic.ttl.variant $DOCKER_DIR
cp -p templ_civic.ttl.drug $DOCKER_DIR
cp -p templ_civic.ttl.disease $DOCKER_DIR
cp -p so_table.tsv $DOCKER_DIR
cp -p drug_table.tsv $DOCKER_DIR
cp -p disease_table.tsv $DOCKER_DIR
sudo docker build -t med2rdf/civic $DOCKER_DIR

