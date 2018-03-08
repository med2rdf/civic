#!/bin/bash
sudo docker run -v `pwd`:/mnt med2rdf/civic -c /mnt/tsv2rdf_civic.json
