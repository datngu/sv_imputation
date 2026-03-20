#!/bin/bash
cd data/array_plink

echo "prefix" > ../input_manifest.csv

for f in *.bed; do
    prefix="${f%.bed}"
    echo "data/array_plink/$prefix"
done >> ../input_manifest.csv