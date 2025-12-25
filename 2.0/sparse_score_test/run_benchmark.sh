#!/bin/bash
# run_benchmark.sh

PLINK2=./bin/plink2
SNP_COUNT=1000000
BENCH_DIR=./benchmark_1m
rm -rf $BENCH_DIR
mkdir -p $BENCH_DIR

# 1. Prepare dataset
echo "Extracting $SNP_COUNT SNPs from g1000_eur..."
$PLINK2 \
  --bfile /n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/g1000_eur \
  --thin-count $SNP_COUNT \
  --make-pgen \
  --out $BENCH_DIR/test \
  --silent

# 2. Generate dense score file (3 traits, 50% sparsity)
echo "Generating synthetic dense score file ($SNP_COUNT SNPs)..."
python3 -c "
import random
import os
import sys

pvar_file = '$BENCH_DIR/test.pvar'
score_file = '$BENCH_DIR/dense.score'

with open(pvar_file, 'r') as f:
    lines = f.readlines()
variants = [l.strip().split() for l in lines if not l.startswith('#')]

with open(score_file, 'w') as f:
    f.write('ID Allele S1 S2 S3\n')
    for v in variants:
        vid = v[2]
        alt = v[4]
        # ~50% sparsity for S1, S2, 0% for S3
        s1 = random.uniform(-1, 1) if random.random() > 0.5 else 0
        s2 = random.uniform(-1, 1) if random.random() > 0.5 else 0
        s3 = random.uniform(-1, 1)
        f.write(f'{vid} {alt} {s1:.4f} {s2:.4f} {s3:.4f}\n')
"

# 3. Convert to sparse
echo "Converting to sparse format..."
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
python3 "$SCRIPT_DIR/convert_to_sparse.py" $BENCH_DIR/dense.score $BENCH_DIR/sparse

# 4. Run benchmarks
echo "Running dense score benchmark..."
/usr/bin/time -v $PLINK2 \
  --pfile $BENCH_DIR/test \
  --score $BENCH_DIR/dense.score 1 2 header-read cols=+scoresums \
  --score-col-nums 3-5 \
  --out $BENCH_DIR/out_dense \
  --threads 2 2>&1 | tee $BENCH_DIR/dense.bench

echo "Running sparse score benchmark..."
/usr/bin/time -v $PLINK2 \
  --pfile $BENCH_DIR/test \
  --score-sparse $BENCH_DIR/sparse.mtx $BENCH_DIR/sparse.snp $BENCH_DIR/sparse.names cols=+scoresums \
  --out $BENCH_DIR/out_sparse \
  --threads 2 2>&1 | tee $BENCH_DIR/sparse.bench

# 5. Verify results
echo "Verifying numerical results..."
python3 -c "
import pandas as pd
import numpy as np
import sys

try:
    dense = pd.read_csv('$BENCH_DIR/out_dense.sscore', sep='\t')
    sparse = pd.read_csv('$BENCH_DIR/out_sparse.sscore', sep='\t')
    
    cols_to_compare = ['NAMED_ALLELE_DOSAGE_SUM', 'S1_SUM', 'S2_SUM', 'S3_SUM']
    for col in cols_to_compare:
        diff = np.abs(dense[col] - sparse[col]).max()
        print(f'Max difference for {col}: {diff}')
        if diff > 1e-2:
            print(f'FAILURE: Numerical mismatch in {col}!')
            sys.exit(1)
    print('SUCCESS: Results match.')
except Exception as e:
    print(f'Error comparing files: {e}')
    sys.exit(1)
"

echo "Benchmark complete. Results in $BENCH_DIR"
