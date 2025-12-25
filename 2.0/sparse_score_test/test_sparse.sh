#!/bin/bash
# test_sparse.sh

PLINK2=./bin/plink2
TEST_DIR=./test_data
rm -rf $TEST_DIR
mkdir -p $TEST_DIR

# 1. Create dummy pgen/pvar/psam
# Extract small subset from g1000_eur (Chromosome 22 for speed)
# Using existing plink2 first to create data? No, use the newly compiled one if it works for basic file ops.
if [ ! -f "$PLINK2" ]; then
    echo "Error: plink2 executable not found at $PLINK2"
    exit 1
fi

echo "Generating test genotype data from g1000_eur..."
$PLINK2 \
  --bfile /n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/g1000_eur \
  --chr 22 \
  --from-bp 16050000 --to-bp 18000000 \
  --make-pgen \
  --out $TEST_DIR/test \
  --silent

# 2. Create dense score file
# We will use python to generate valid scores matching the PVAR
echo "Generating synthetic dense score file..."
python3 -c "
import random
with open('$TEST_DIR/test.pvar', 'r') as f:
    lines = f.readlines()
    
# Skip header lines (start with #)
variants = [l.strip().split() for l in lines if not l.startswith('#')]

with open('$TEST_DIR/dense.score', 'w') as f:
    f.write('ID Allele S1 S2 S3\n')
    for v in variants:
        # ID is col 2 (0-indexed), Ref col 3, Alt col 4.
        # PLINK pvar: #CHROM POS ID REF ALT
        # We'll use ALT allele for scoring (col 4)
        vid = v[2]
        alt = v[4]
        # Random weights, sometimes 0 to test sparsity
        s1 = random.uniform(-1, 1) if random.random() > 0.5 else 0
        s2 = random.uniform(-1, 1) if random.random() > 0.8 else 0
        s3 = random.uniform(-1, 1)
        f.write(f'{vid} {alt} {s1:.4f} {s2:.4f} {s3:.4f}\n')
"

# 3. Convert to sparse
echo "Converting dense score to sparse format..."
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python3 "$SCRIPT_DIR/convert_to_sparse.py" $TEST_DIR/dense.score $TEST_DIR/sparse

# 4. Run dense score
echo "Running dense score calculation..."
$PLINK2 \
  --pfile $TEST_DIR/test \
  --score $TEST_DIR/dense.score 1 2 header-read cols=+scoresums \
  --out $TEST_DIR/out_dense \
  --threads 2

# 5. Run sparse score
echo "Running sparse score calculation..."
$PLINK2 \
  --pfile $TEST_DIR/test \
  --score-sparse $TEST_DIR/sparse.mtx $TEST_DIR/sparse.snp $TEST_DIR/sparse.names cols=+scoresums \
  --out $TEST_DIR/out_sparse \
  --threads 2

# 6. Compare
echo "Comparing results..."
# Use python to compare numerical values
python3 -c "
import pandas as pd
import numpy as np
import sys

try:
    dense = pd.read_csv('$TEST_DIR/out_dense.sscore', sep='\t')
    sparse = pd.read_csv('$TEST_DIR/out_sparse.sscore', sep='\t')
    
    # Check if headers match (sparse might name differently?)
    # Dense: ID, IID, ... S1_SUM, S2_SUM, S3_SUM
    # Sparse: My implementation outputs score names as given.
    # If dense uses header-read, it uses names from file.
    
    # Columns to compare: S1_SUM, S2_SUM, S3_SUM
    cols = ['S1_SUM', 'S2_SUM', 'S3_SUM']
    
    # Note: Column names in standard plink output might append _SUM if requested.
    # My sparse implementation might default to _SUM or AVG?
    # I set 'cols=+scoresums' for dense.
    # For sparse, I assume my implementation defaults to something.
    # Let's check output columns.
    
    # Check common columns
    common_cols = [c for c in dense.columns if c in sparse.columns and 'SUM' in c]
    if not common_cols:
        print('No common score columns found!')
        print('Dense cols:', dense.columns)
        print('Sparse cols:', sparse.columns)
        sys.exit(1)
        
    print(f'Comparing columns: {common_cols}')
    
    diff_max = 0
    for c in common_cols:
        d = dense[c].values
        s = sparse[c].values
        diff = np.abs(d - s)
        max_d = np.max(diff)
        if max_d > diff_max: diff_max = max_d
        
    print(f'Max difference: {diff_max}')
    if diff_max < 1e-4:
        print('SUCCESS: Outputs match within tolerance.')
    else:
        print('FAILURE: Outputs differ.')
        sys.exit(1)
        
except Exception as e:
    print(f'Error comparing files: {e}')
    sys.exit(1)
"
