#!/bin/bash

OUT_DIR='output'
MAF=example_maf.txt
BIN_MATRIX="$OUT_DIR"/bin_mat.txt
INT_MATRIX="$OUT_DIR"/int_mat.txt
mkdir -p $OUT_DIR

python process_maf.py -maf $MAF -bo $BIN_MATRIX -io $INT_MATRIX