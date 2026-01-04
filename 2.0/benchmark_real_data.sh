#!/bin/bash
# Benchmark Sparse Scoring vs Traditional PLINK2 Scoring
# Using real genotype data from 1000 Genomes EUR
# 1M SNPs × 100 scores with 50% sparsity

set -e

# Configuration
GENO_DATA="/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/g1000_eur"
BENCH_DIR="./benchmark_1m_real"
NUM_SNPS=1000000
NUM_SCORES=100
SPARSITY=0.50

# Try to find plink2
if [ -f "./bin/plink2" ]; then
    PLINK2="./bin/plink2"
elif [ -f "/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/plink2" ]; then
    PLINK2="/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/plink2"
elif command -v plink2 &> /dev/null; then
    PLINK2="plink2"
elif [ -f "./plink2" ]; then
    PLINK2="./plink2"
else
    echo "Error: plink2 not found"
    echo "Please build plink2 or specify PLINK2 path"
    exit 1
fi

echo "=========================================="
echo "  Real Data Benchmark: Sparse vs Dense"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Genotype data: ${GENO_DATA}"
echo "  SNPs: ${NUM_SNPS}"
echo "  Scores: ${NUM_SCORES}"
echo "  Sparsity: ${SPARSITY}"
echo "  PLINK2: ${PLINK2}"
echo ""

# Create benchmark directory
mkdir -p $BENCH_DIR

# Step 1: Extract 1M SNPs from real data
echo "=========================================="
echo "  Step 1: Extracting ${NUM_SNPS} SNPs"
echo "=========================================="

if [ ! -f "${BENCH_DIR}/data.pgen" ]; then
    echo "Extracting SNPs from ${GENO_DATA}..."
    $PLINK2 --bfile $GENO_DATA \
        --thin-count $NUM_SNPS \
        --make-pgen \
        --out ${BENCH_DIR}/data \
        --threads 4 \
        2>&1 | tee ${BENCH_DIR}/extract.log
    
    echo "✓ Genotype data extracted"
else
    echo "✓ Genotype data already exists"
fi

# Get actual number of variants extracted
ACTUAL_SNPS=$(wc -l < ${BENCH_DIR}/data.pvar | tail -1)
if [ -f "${BENCH_DIR}/data.pvar" ]; then
    ACTUAL_SNPS=$(grep -v "^#" ${BENCH_DIR}/data.pvar | wc -l)
else
    ACTUAL_SNPS=$NUM_SNPS
fi
echo "  Actual SNPs extracted: ${ACTUAL_SNPS}"

# Step 2: Generate sparse score matrix
echo ""
echo "=========================================="
echo "  Step 2: Generating Sparse Score Matrix"
echo "=========================================="

python3 - <<EOF
import random
import numpy as np
import struct
import sys

num_variants = ${ACTUAL_SNPS}
num_traits = ${NUM_SCORES}
sparsity = ${SPARSITY}
bench_dir = "${BENCH_DIR}"

print(f"Generating sparse matrix:")
print(f"  Variants: {num_variants}")
print(f"  Traits: {num_traits}")
print(f"  Sparsity: {sparsity:.1%}")

# Read variant IDs from pvar file
variant_ids = []
try:
    with open(f'{bench_dir}/data.pvar', 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                variant_ids.append(parts[2])
except Exception as e:
    print(f"Error reading pvar: {e}")
    sys.exit(1)

num_variants = min(num_variants, len(variant_ids))
variant_ids = variant_ids[:num_variants]

print(f"Using {num_variants} variants")

# Generate sparse matrix
np.random.seed(42)
row_starts = [0]
col_indices = []
values = []
nonzeros = 0

print("Generating sparse weights...")
for variant_idx in range(num_variants):
    if variant_idx > 0 and variant_idx % 100000 == 0:
        print(f"  Processed {variant_idx}/{num_variants} variants...")
    
    for trait_idx in range(num_traits):
        if random.random() > sparsity:
            col_indices.append(trait_idx)
            values.append(np.random.normal(0, 1))
            nonzeros += 1
    row_starts.append(nonzeros)

density = nonzeros / (num_variants * num_traits)
print(f"Generated {nonzeros} non-zero weights")
print(f"Actual density: {density:.2%}")

# Write MTX format
print("Writing MTX format...")
with open(f"{bench_dir}/sparse.mtx", 'w') as f:
    f.write("%%MatrixMarket matrix coordinate real general\n")
    f.write(f"{num_variants} {num_traits} {nonzeros}\n")
    for variant_idx in range(num_variants):
        start = row_starts[variant_idx]
        end = row_starts[variant_idx + 1]
        for k in range(start, end):
            f.write(f"{variant_idx + 1} {col_indices[k] + 1} {values[k]:.6f}\n")

# Write SNP file
print("Writing SNP file...")
with open(f"{bench_dir}/sparse.snp", 'w') as f:
    for vid in variant_ids:
        f.write(f"{vid}\n")

# Write names file
print("Writing names file...")
with open(f"{bench_dir}/sparse.names", 'w') as f:
    for i in range(num_traits):
        f.write(f"Score{i+1}\n")

# Write binary format
print("Writing binary format...")
with open(f"{bench_dir}/sparse.binsparse", 'wb') as f:
    # Header
    f.write(b'PLNKSPR\x00')  # Magic
    f.write(struct.pack('<I', 1))  # Version
    f.write(struct.pack('<QQQ', num_variants, num_traits, nonzeros))  # Dimensions
    f.write(b'\x00' * 32)  # Reserved
    
    # CSR arrays
    row_starts_array = np.array(row_starts, dtype=np.uint64)
    f.write(row_starts_array.tobytes())
    
    col_indices_array = np.array(col_indices, dtype=np.uint32)
    f.write(col_indices_array.tobytes())
    
    values_array = np.array(values, dtype=np.float64)
    f.write(values_array.tobytes())
    
    # Variant IDs
    for vid in variant_ids:
        f.write(vid.encode('utf-8') + b'\x00')
    
    # Trait names
    for i in range(num_traits):
        f.write(f"Score{i+1}".encode('utf-8') + b'\x00')

print(f"✓ Sparse matrix files created")
print(f"  - {bench_dir}/sparse.mtx")
print(f"  - {bench_dir}/sparse.binsparse")
print(f"  - {bench_dir}/sparse.snp")
print(f"  - {bench_dir}/sparse.names")
EOF

# Step 3: Generate dense score file for comparison
echo ""
echo "=========================================="
echo "  Step 3: Generating Dense Score File"
echo "=========================================="

python3 - <<EOF
import random
import numpy as np
import sys

num_variants = ${ACTUAL_SNPS}
num_traits = ${NUM_SCORES}
sparsity = ${SPARSITY}
bench_dir = "${BENCH_DIR}"

print(f"Generating dense score file...")

# Read variant IDs
variant_ids = []
with open(f'{bench_dir}/data.pvar', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 3:
            variant_ids.append(parts[2])

num_variants = min(num_variants, len(variant_ids))
variant_ids = variant_ids[:num_variants]

# Read sparse matrix to create equivalent dense
row_starts = []
col_indices = []
values = []

with open(f'{bench_dir}/sparse.mtx', 'r') as f:
    # Skip header
    f.readline()
    f.readline()
    
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3:
            variant_idx = int(parts[0]) - 1
            trait_idx = int(parts[1]) - 1
            value = float(parts[2])
            
            # Expand row_starts if needed
            while len(row_starts) <= variant_idx + 1:
                row_starts.append(len(col_indices))
            
            col_indices.append(trait_idx)
            values.append(value)

# Fill in final row_starts
while len(row_starts) <= num_variants:
    row_starts.append(len(col_indices))

# Create dense matrix
dense_matrix = np.zeros((num_variants, num_traits))
for variant_idx in range(num_variants):
    if variant_idx > 0 and variant_idx % 100000 == 0:
        print(f"  Processed {variant_idx}/{num_variants} variants...")
    
    start = row_starts[variant_idx]
    end = row_starts[min(variant_idx + 1, len(row_starts) - 1)]
    
    for k in range(start, end):
        trait = col_indices[k]
        weight = values[k]
        dense_matrix[variant_idx, trait] = weight

# Write dense score file
print("Writing dense score file...")
with open(f'{bench_dir}/dense.score', 'w') as f:
    # Header
    f.write('SNP\tA1')
    for i in range(num_traits):
        f.write(f'\tS{i+1}')
    f.write('\n')
    
    # Data
    for idx, vid in enumerate(variant_ids):
        if idx > 0 and idx % 100000 == 0:
            print(f"  Writing {idx}/{num_variants} variants...")
        
        f.write(f'{vid}\tA')
        for trait_idx in range(num_traits):
            f.write(f'\t{dense_matrix[idx, trait_idx]:.6f}')
        f.write('\n')

print(f"✓ Dense score file created: {bench_dir}/dense.score")
EOF

# Step 4: Run benchmarks
echo ""
echo "=========================================="
echo "  Step 4: Running Benchmarks"
echo "=========================================="

# Benchmark 1: Traditional Dense Scoring
echo ""
echo "--- Benchmark 1: Traditional Dense Scoring ---"
echo ""

DENSE_START=$(date +%s.%N)
$PLINK2 --pfile ${BENCH_DIR}/data \
    --score ${BENCH_DIR}/dense.score 1 2 header-read cols=+scoresums \
    --score-col-nums 3-$((2+${NUM_SCORES})) \
    --out ${BENCH_DIR}/results_dense \
    --threads 4 \
    2>&1 | tee ${BENCH_DIR}/benchmark_dense.log
DENSE_END=$(date +%s.%N)
DENSE_TIME=$(echo "$DENSE_END - $DENSE_START" | bc)

echo ""
echo "Dense scoring completed in ${DENSE_TIME}s"

# Check if sparse scoring is supported
# Note: --score-sparse may work even if not in help text, so test it directly
SPARSE_SUPPORTED=false
SPARSE_BIN_SUCCESS=false
SPARSE_MTX_SUCCESS=false
SPARSE_BIN_TIME=0
SPARSE_MTX_TIME=0

# First check if we have the required sparse files
if [ ! -f "${BENCH_DIR}/sparse.mtx" ] || [ ! -f "${BENCH_DIR}/sparse.snp" ] || [ ! -f "${BENCH_DIR}/sparse.names" ]; then
    echo "⚠ Sparse matrix files not found"
    echo "  Skipping sparse scoring benchmarks"
    SPARSE_SUPPORTED=false
else
    # Test if --score-sparse is recognized by trying to parse it
    # If it's unrecognized, plink2 will error immediately with "Unrecognized flag"
    TEST_OUTPUT=$($PLINK2 --pfile ${BENCH_DIR}/data --score-sparse ${BENCH_DIR}/sparse.mtx ${BENCH_DIR}/sparse.snp ${BENCH_DIR}/sparse.names --out ${BENCH_DIR}/sparse_test_check 2>&1)
    
    if echo "$TEST_OUTPUT" | grep -qE "(Error.*Unrecognized.*score-sparse|Error.*score-sparse.*not found)"; then
        SPARSE_SUPPORTED=false
        echo "⚠ Sparse scoring is NOT supported in this PLINK2 version"
        echo "  Skipping sparse scoring benchmarks"
        echo "  Note: Sparse scoring requires a custom PLINK2 build with sparse scoring support"
        # Clean up test files
        rm -f ${BENCH_DIR}/sparse_test_check.* 2>/dev/null
    else
        SPARSE_SUPPORTED=true
        echo "✓ Sparse scoring is supported in this PLINK2 version"
        # Keep test files in benchmark directory (don't delete)
    fi
fi

# Benchmark 2: Sparse Scoring (Binary Format)
if [ "$SPARSE_SUPPORTED" = true ]; then
    echo ""
    echo "--- Benchmark 2: Sparse Scoring (Binary Format) ---"
    echo ""

    SPARSE_BIN_START=$(date +%s.%N)
    if $PLINK2 --pfile ${BENCH_DIR}/data \
        --score-sparse ${BENCH_DIR}/sparse.binsparse ${BENCH_DIR}/sparse.snp ${BENCH_DIR}/sparse.names \
        cols=+scoresums \
        --out ${BENCH_DIR}/results_sparse_binary \
        --threads 4 \
        2>&1 | tee ${BENCH_DIR}/benchmark_sparse_binary.log; then
        SPARSE_BIN_END=$(date +%s.%N)
        SPARSE_BIN_TIME=$(echo "$SPARSE_BIN_END - $SPARSE_BIN_START" | bc)
        echo ""
        echo "Sparse binary scoring completed in ${SPARSE_BIN_TIME}s"
        SPARSE_BIN_SUCCESS=true
    else
        echo ""
        echo "✗ Sparse binary scoring failed"
        SPARSE_BIN_TIME=0
        SPARSE_BIN_SUCCESS=false
    fi

    # Benchmark 3: Sparse Scoring (MTX Format)
    echo ""
    echo "--- Benchmark 3: Sparse Scoring (MTX Format) ---"
    echo ""

    SPARSE_MTX_START=$(date +%s.%N)
    if $PLINK2 --pfile ${BENCH_DIR}/data \
        --score-sparse ${BENCH_DIR}/sparse.mtx ${BENCH_DIR}/sparse.snp ${BENCH_DIR}/sparse.names \
        cols=+scoresums \
        --out ${BENCH_DIR}/results_sparse_mtx \
        --threads 4 \
        2>&1 | tee ${BENCH_DIR}/benchmark_sparse_mtx.log; then
        SPARSE_MTX_END=$(date +%s.%N)
        SPARSE_MTX_TIME=$(echo "$SPARSE_MTX_END - $SPARSE_MTX_START" | bc)
        echo ""
        echo "Sparse MTX scoring completed in ${SPARSE_MTX_TIME}s"
        SPARSE_MTX_SUCCESS=true
    else
        echo ""
        echo "✗ Sparse MTX scoring failed"
        SPARSE_MTX_TIME=0
        SPARSE_MTX_SUCCESS=false
    fi
fi

# Step 5: Compare results
echo ""
echo "=========================================="
echo "  Step 5: Comparing Results"
echo "=========================================="

python3 - <<EOF
import pandas as pd
import numpy as np
import sys
import os

bench_dir = "${BENCH_DIR}"
sparse_supported = "${SPARSE_SUPPORTED}" == "true"
sparse_bin_success = "${SPARSE_BIN_SUCCESS:-false}" == "true"
sparse_mtx_success = "${SPARSE_MTX_SUCCESS:-false}" == "true"

print("Comparing scoring results...")

# Read dense score file
try:
    if os.path.exists(f'{bench_dir}/results_dense.sscore'):
        df_dense = pd.read_csv(f'{bench_dir}/results_dense.sscore', sep='\t', comment='#')
        score_cols_dense = [c for c in df_dense.columns if c.startswith('SCORE')]
        print(f"\nDense scoring: {len(score_cols_dense)} score columns found")
    else:
        print(f"\n✗ Dense score file not found: {bench_dir}/results_dense.sscore")
        df_dense = None
        score_cols_dense = []
except Exception as e:
    print(f"Error reading dense scores: {e}")
    df_dense = None
    score_cols_dense = []

# Compare dense vs sparse binary
if sparse_bin_success and df_dense is not None:
    try:
        if os.path.exists(f'{bench_dir}/results_sparse_binary.sscore'):
            df_sparse_bin = pd.read_csv(f'{bench_dir}/results_sparse_binary.sscore', sep='\t', comment='#')
            score_cols_bin = [c for c in df_sparse_bin.columns if c.startswith('SCORE')]
            
            print(f"\nSparse Binary: {len(score_cols_bin)} score columns found")
            
            if len(score_cols_dense) == len(score_cols_bin):
                max_diff_bin = 0.0
                for col_d, col_b in zip(score_cols_dense, score_cols_bin):
                    diff = np.abs(df_dense[col_d] - df_sparse_bin[col_b])
                    max_diff_bin = max(max_diff_bin, diff.max())
                
                print(f"Dense vs Sparse Binary:")
                print(f"  Max difference: {max_diff_bin:.2e}")
                if max_diff_bin < 1e-5:
                    print(f"  ✓ Results match (within numerical precision)")
                else:
                    print(f"  ✗ Results differ")
            else:
                print(f"  ⚠ Column count mismatch: {len(score_cols_dense)} vs {len(score_cols_bin)}")
        else:
            print(f"\n✗ Sparse binary score file not found")
    except Exception as e:
        print(f"Error comparing sparse binary: {e}")
elif not sparse_supported:
    print(f"\n⚠ Sparse binary comparison skipped (not supported)")
else:
    print(f"\n⚠ Sparse binary comparison skipped (scoring failed)")

# Compare dense vs sparse MTX
if sparse_mtx_success and df_dense is not None:
    try:
        if os.path.exists(f'{bench_dir}/results_sparse_mtx.sscore'):
            df_sparse_mtx = pd.read_csv(f'{bench_dir}/results_sparse_mtx.sscore', sep='\t', comment='#')
            score_cols_mtx = [c for c in df_sparse_mtx.columns if c.startswith('SCORE')]
            
            print(f"\nSparse MTX: {len(score_cols_mtx)} score columns found")
            
            if len(score_cols_dense) == len(score_cols_mtx):
                max_diff_mtx = 0.0
                for col_d, col_m in zip(score_cols_dense, score_cols_mtx):
                    diff = np.abs(df_dense[col_d] - df_sparse_mtx[col_m])
                    max_diff_mtx = max(max_diff_mtx, diff.max())
                
                print(f"Dense vs Sparse MTX:")
                print(f"  Max difference: {max_diff_mtx:.2e}")
                if max_diff_mtx < 1e-5:
                    print(f"  ✓ Results match (within numerical precision)")
                else:
                    print(f"  ✗ Results differ")
            else:
                print(f"  ⚠ Column count mismatch: {len(score_cols_dense)} vs {len(score_cols_mtx)}")
        else:
            print(f"\n✗ Sparse MTX score file not found")
    except Exception as e:
        print(f"Error comparing sparse MTX: {e}")
elif not sparse_supported:
    print(f"\n⚠ Sparse MTX comparison skipped (not supported)")
else:
    print(f"\n⚠ Sparse MTX comparison skipped (scoring failed)")
EOF

# Step 6: Performance Summary
echo ""
echo "=========================================="
echo "  Performance Summary"
echo "=========================================="

python3 - <<EOF
import os

dense_time = ${DENSE_TIME}
sparse_bin_time = ${SPARSE_BIN_TIME:-0}
sparse_mtx_time = ${SPARSE_MTX_TIME:-0}
sparse_supported = "${SPARSE_SUPPORTED}" == "true"
sparse_bin_success = "${SPARSE_BIN_SUCCESS:-false}" == "true"
sparse_mtx_success = "${SPARSE_MTX_SUCCESS:-false}" == "true"

print(f"\n{'Method':<25} {'Time (s)':<15} {'Speedup':<15} {'Status':<15}")
print("-" * 70)
print(f"{'Dense (Traditional)':<25} {dense_time:<15.2f} {'1.00x':<15} {'✓':<15}")

if sparse_supported and sparse_bin_success and sparse_bin_time > 0:
    speedup = dense_time / sparse_bin_time if sparse_bin_time > 0 else 0
    print(f"{'Sparse (Binary)':<25} {sparse_bin_time:<15.2f} {speedup:<15.2f}x {'✓':<15}")
elif sparse_supported:
    print(f"{'Sparse (Binary)':<25} {'N/A':<15} {'N/A':<15} {'✗ Failed':<15}")
else:
    print(f"{'Sparse (Binary)':<25} {'N/A':<15} {'N/A':<15} {'Not supported':<15}")

if sparse_supported and sparse_mtx_success and sparse_mtx_time > 0:
    speedup = dense_time / sparse_mtx_time if sparse_mtx_time > 0 else 0
    print(f"{'Sparse (MTX)':<25} {sparse_mtx_time:<15.2f} {speedup:<15.2f}x {'✓':<15}")
elif sparse_supported:
    print(f"{'Sparse (MTX)':<25} {'N/A':<15} {'N/A':<15} {'✗ Failed':<15}")
else:
    print(f"{'Sparse (MTX)':<25} {'N/A':<15} {'N/A':<15} {'Not supported':<15}")

print(f"\nMemory efficiency:")
bench_dir = "${BENCH_DIR}"
actual_snps = ${ACTUAL_SNPS}
num_scores = ${NUM_SCORES}

dense_size_gb = (actual_snps * num_scores * 8) / (1024**3)
print(f"  Dense format: ~{dense_size_gb:.2f} GB")

if os.path.exists(f'{bench_dir}/sparse.binsparse'):
    sparse_size = os.path.getsize(f'{bench_dir}/sparse.binsparse')
    sparse_size_gb = sparse_size / (1024**3)
    print(f"  Sparse format: ~{sparse_size_gb:.2f} GB")
    print(f"  Memory savings: {(1 - sparse_size_gb/dense_size_gb)*100:.1f}%")
else:
    print(f"  Sparse format: File not found")
EOF

echo ""
echo "=========================================="
echo "  Benchmark Complete"
echo "=========================================="
echo ""
echo "Results saved in: ${BENCH_DIR}/"
echo "  - results_dense.sscore"

if [ "$SPARSE_SUPPORTED" = true ]; then
    if [ "$SPARSE_BIN_SUCCESS" = true ]; then
        echo "  - results_sparse_binary.sscore"
    else
        echo "  - results_sparse_binary.sscore (failed)"
    fi
    
    if [ "$SPARSE_MTX_SUCCESS" = true ]; then
        echo "  - results_sparse_mtx.sscore"
    else
        echo "  - results_sparse_mtx.sscore (failed)"
    fi
else
    echo "  - Sparse scoring not available in this PLINK2 version"
    echo ""
    echo "Note: To test sparse scoring, you need a PLINK2 build with sparse scoring support."
    echo "      The sparse score matrices have been generated and are available at:"
    echo "      - ${BENCH_DIR}/sparse.mtx"
    echo "      - ${BENCH_DIR}/sparse.binsparse"
    echo "      - ${BENCH_DIR}/sparse.snp"
    echo "      - ${BENCH_DIR}/sparse.names"
fi
echo ""

