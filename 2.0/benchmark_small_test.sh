#!/bin/bash
# Small Benchmark Test: Sparse Scoring vs Traditional PLINK2 Scoring
# Using real genotype data from 1000 Genomes EUR
# 1000 SNPs × 100 scores with 50% sparsity

set -e

# Configuration
GENO_DATA="/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/g1000_eur"
NUM_SNPS=6000000
NUM_SCORES=100
# 1M non-zero per score out of 6M = 1/6 density = 83.3% sparsity
SPARSITY=0.8333
BENCH_DIR="./benchmark_${NUM_SNPS}_snps"

# Use system plink2 for dense scoring (baseline/original version)
if [ -f "/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/plink2" ]; then
    PLINK2_SYSTEM="/n/holylfs05/LABS/liang_lab/Lab/btruong/Tools/plink2"
elif command -v plink2 &> /dev/null; then
    PLINK2_SYSTEM="plink2"
else
    echo "Error: System plink2 not found"
    exit 1
fi

# Use built plink2 for sparse scoring (with sparse scoring support)
if [ -f "./bin/plink2" ]; then
    PLINK2_BUILT="./bin/plink2"
    echo "✓ Found built plink2 at ./bin/plink2"
else
    echo "⚠ Built plink2 not found at ./bin/plink2"
    echo "  Please build plink2 first with: make plink2"
    echo "  Or run: bash build_plink2.sh"
    PLINK2_BUILT=""
fi

# Default to system plink2 for dense scoring
PLINK2=$PLINK2_SYSTEM

echo "=========================================="
echo "  Small Test: Sparse vs Dense Scoring"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Genotype data: ${GENO_DATA}"
echo "  SNPs: ${NUM_SNPS}"
echo "  Scores: ${NUM_SCORES}"
echo "  Sparsity: ${SPARSITY}"
echo "  System PLINK2 (dense): ${PLINK2_SYSTEM}"
if [ -n "$PLINK2_BUILT" ]; then
    echo "  Built PLINK2 (sparse): ${PLINK2_BUILT}"
else
    echo "  Built PLINK2 (sparse): Not available"
fi
echo ""

# Create benchmark directory
mkdir -p $BENCH_DIR

# Step 1: Extract 1K SNPs from real data
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
    if variant_idx > 0 and variant_idx % 200 == 0:
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
    # Header (must match C struct alignment: 4 bytes padding after version)
    f.write(b'PLNKSPR\x00')  # Magic (8 bytes)
    f.write(struct.pack('<I', 1))  # Version (4 bytes)
    f.write(b'\x00' * 4)  # Padding for uint64_t alignment (4 bytes)
    f.write(struct.pack('<QQQ', num_variants, num_traits, nonzeros))  # Dimensions (24 bytes)
    f.write(b'\x00' * 32)  # Reserved (32 bytes)
    
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
# First, build a dictionary mapping variant_idx -> {trait_idx: value}
sparse_dict = {}
size_line_seen = False
with open(f'{bench_dir}/sparse.mtx', 'r') as f:
    # Read all lines and parse
    for line in f:
        line = line.strip()
        # Skip comments
        if line.startswith('%'):
            continue
        parts = line.split()
        if len(parts) >= 3:
            try:
                r = int(parts[0])
                c = int(parts[1])
                v = float(parts[2])
                # Skip size line (first non-comment line)
                if not size_line_seen:
                    size_line_seen = True
                    continue
                variant_idx = r - 1  # Convert 1-based to 0-based
                trait_idx = c - 1    # Convert 1-based to 0-based
                if variant_idx not in sparse_dict:
                    sparse_dict[variant_idx] = {}
                sparse_dict[variant_idx][trait_idx] = v
            except ValueError:
                continue

# Create dense matrix from sparse dictionary
dense_matrix = np.zeros((num_variants, num_traits))
for variant_idx in range(num_variants):
    if variant_idx in sparse_dict:
        for trait_idx, value in sparse_dict[variant_idx].items():
            if trait_idx < num_traits:
                dense_matrix[variant_idx, trait_idx] = value

# Read alleles from pvar file
alleles = {}
with open(f'{bench_dir}/data.pvar', 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            vid = parts[2]
            alt_allele = parts[4]  # Use ALT allele (column 5, 0-indexed as 4)
            alleles[vid] = alt_allele

# Write dense score file with correct alleles
print("Writing dense score file...")
with open(f'{bench_dir}/dense.score', 'w') as f:
    # Header
    f.write('SNP\tA1')
    for i in range(num_traits):
        f.write(f'\tS{i+1}')
    f.write('\n')
    
    # Data
    for idx, vid in enumerate(variant_ids):
        allele = alleles.get(vid, 'A')
        f.write(f'{vid}\t{allele}')
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

# Initialize variables
SPARSE_SUPPORTED=false
SPARSE_BIN_SUCCESS=false
SPARSE_MTX_SUCCESS=false
SPARSE_BIN_TIME=0
SPARSE_MTX_TIME=0

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

# Check if sparse scoring is supported in built plink2
# Note: --score-sparse works even if not in help text, so we test it directly
if [ -n "$PLINK2_BUILT" ] && [ -f "$PLINK2_BUILT" ]; then
    # Test if --score-sparse is recognized by trying to parse it
    # If it's unrecognized, plink2 will error immediately with "Unrecognized flag"
    TEST_OUTPUT=$($PLINK2_BUILT --pfile ${BENCH_DIR}/data --score-sparse ${BENCH_DIR}/sparse.mtx ${BENCH_DIR}/sparse.snp ${BENCH_DIR}/sparse.names --out ${BENCH_DIR}/sparse_test_check 2>&1)
    
    if echo "$TEST_OUTPUT" | grep -qE "(Error.*Unrecognized.*score-sparse|Error.*score-sparse.*not found)"; then
        SPARSE_SUPPORTED=false
        echo "⚠ Sparse scoring is NOT supported in built PLINK2 version"
        echo "  Skipping sparse scoring benchmarks"
        # Clean up test files
        rm -f ${BENCH_DIR}/sparse_test_check.* 2>/dev/null
    else
        SPARSE_SUPPORTED=true
        echo "✓ Sparse scoring is supported in built PLINK2 version"
        # Keep test files in benchmark directory (don't delete)
    fi
else
    echo "⚠ Built PLINK2 not available"
    echo "  Skipping sparse scoring benchmarks"
    echo "  Note: Please build plink2 with sparse scoring support"
    SPARSE_SUPPORTED=false
fi

# Benchmark 2: Sparse Scoring (Binary Format) - Use built plink2
if [ "$SPARSE_SUPPORTED" = true ]; then
    echo ""
    echo "--- Benchmark 2: Sparse Scoring (Binary Format) ---"
    echo "  Using built PLINK2: ${PLINK2_BUILT}"
    echo ""

    SPARSE_BIN_START=$(date +%s.%N)
    if $PLINK2_BUILT --pfile ${BENCH_DIR}/data \
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

    # Benchmark 3: Sparse Scoring (MTX Format) - Use built plink2
    echo ""
    echo "--- Benchmark 3: Sparse Scoring (MTX Format) ---"
    echo "  Using built PLINK2: ${PLINK2_BUILT}"
    echo ""

    SPARSE_MTX_START=$(date +%s.%N)
    if $PLINK2_BUILT --pfile ${BENCH_DIR}/data \
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
import sys
import os

# Try to import pandas, but handle import errors gracefully
try:
    import pandas as pd
    import numpy as np
    PANDAS_AVAILABLE = True
except ImportError as e:
    print(f"Warning: pandas not available ({e})")
    print("Skipping detailed comparison - will show file existence only")
    PANDAS_AVAILABLE = False

bench_dir = "${BENCH_DIR}"
sparse_supported = "${SPARSE_SUPPORTED}" == "true"
sparse_bin_success = "${SPARSE_BIN_SUCCESS:-false}" == "true"
sparse_mtx_success = "${SPARSE_MTX_SUCCESS:-false}" == "true"

print("Comparing scoring results...")

if not PANDAS_AVAILABLE:
    # Simple file-based comparison
    print("\nFile status:")
    if os.path.exists(f'{bench_dir}/results_dense.sscore'):
        print(f"  ✓ Dense score file exists")
    else:
        print(f"  ✗ Dense score file not found")
    
    if sparse_supported and sparse_bin_success:
        if os.path.exists(f'{bench_dir}/results_sparse_binary.sscore'):
            print(f"  ✓ Sparse binary score file exists")
        else:
            print(f"  ✗ Sparse binary score file not found")
    
    if sparse_supported and sparse_mtx_success:
        if os.path.exists(f'{bench_dir}/results_sparse_mtx.sscore'):
            print(f"  ✓ Sparse MTX score file exists")
        else:
            print(f"  ✗ Sparse MTX score file not found")
    
    sys.exit(0)

# Read dense score file (ground truth from system plink2)
print("\nUsing system PLINK2 dense scoring results as ground truth...")
try:
    if os.path.exists(f'{bench_dir}/results_dense.sscore'):
        # Handle header that starts with #
        with open(f'{bench_dir}/results_dense.sscore', 'r') as f:
            header = f.readline().strip()
            if header.startswith('#'):
                header = header[1:]  # Remove the #
            cols = header.split('\t')
        
        # Read data with proper column names
        df_dense = pd.read_csv(f'{bench_dir}/results_dense.sscore', sep='\t', comment='#', names=cols, skiprows=1)
        
        # Find score columns - can be S1_AVG, S2_AVG, etc. or Score1_AVG, etc.
        # Match pattern: contains number followed by _AVG (exclude NAMED_ALLELE_DOSAGE_SUM)
        import re
        score_cols_dense = [c for c in cols if '_AVG' in c and c != 'NAMED_ALLELE_DOSAGE_SUM' and re.search(r'\d+', c)]
        # Sort by number to ensure consistent ordering
        def extract_num(col):
            match = re.search(r'(\d+)', col)
            return int(match.group(1)) if match else 0
        score_cols_dense = sorted(score_cols_dense, key=extract_num)
        print(f"  Ground truth (dense): {len(score_cols_dense)} score columns found")
        print(f"  Samples: {len(df_dense)}")
        print(f"  Score columns: {score_cols_dense[:5]}..." if len(score_cols_dense) > 5 else f"  Score columns: {score_cols_dense}")
    else:
        print(f"\n✗ Ground truth dense score file not found: {bench_dir}/results_dense.sscore")
        df_dense = None
        score_cols_dense = []
except Exception as e:
    print(f"Error reading ground truth dense scores: {e}")
    import traceback
    traceback.print_exc()
    df_dense = None
    score_cols_dense = []

# Compare dense vs sparse binary (using ground truth)
if sparse_bin_success and df_dense is not None:
    try:
        if os.path.exists(f'{bench_dir}/results_sparse_binary.sscore'):
            # Handle header that starts with # and potentially corrupted/duplicate column names
            with open(f'{bench_dir}/results_sparse_binary.sscore', 'r', encoding='latin-1') as f:
                header = f.readline().strip()
                if header.startswith('#'):
                    header = header[1:]
                cols_bin = header.split('\t')
            
            # Handle duplicate column names by appending index
            seen = {}
            unique_cols = []
            for col in cols_bin:
                if col in seen:
                    seen[col] += 1
                    unique_cols.append(f"{col}_{seen[col]}")
                else:
                    seen[col] = 0
                    unique_cols.append(col)
            cols_bin = unique_cols
            
            df_sparse_bin = pd.read_csv(f'{bench_dir}/results_sparse_binary.sscore', sep='\t', comment='#', names=cols_bin, skiprows=1, encoding='latin-1')
            
            # Find score columns - binary output has corrupted names, so match by position
            # Standard sscore columns: #FID, IID, ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, then AVG cols, then SUM cols
            # So AVG columns start at index 4, and there are len(score_cols_dense) of them
            num_scores = len(score_cols_dense)
            if num_scores > 0 and len(cols_bin) >= 4 + num_scores:
                score_cols_bin = cols_bin[4:4+num_scores]  # Positional match for AVG columns
            else:
                # Fallback: try pattern matching
                import re
                score_cols_bin = [c for c in cols_bin if '_AVG' in c and re.search(r'\d+', c)]
                def extract_num(col):
                    match = re.search(r'(\d+)', col)
                    return int(match.group(1)) if match else 0
                score_cols_bin = sorted(score_cols_bin, key=extract_num)
            
            print(f"\nSparse Binary: {len(score_cols_bin)} score columns found (positional match)")
            
            if len(score_cols_dense) == len(score_cols_bin):
                # Match by index position (ground truth comparison)
                max_diff_bin = 0.0
                mean_diff_bin = 0.0
                total_comparisons = 0
                
                for idx, (col_d, col_b) in enumerate(zip(score_cols_dense, score_cols_bin)):
                    # Match samples by IID
                    if 'IID' in df_dense.columns and 'IID' in df_sparse_bin.columns:
                        merged = pd.merge(df_dense[['IID', col_d]], df_sparse_bin[['IID', col_b]], on='IID', how='inner')
                        if len(merged) > 0:
                            diff = np.abs(merged[col_d] - merged[col_b])
                            max_diff_bin = max(max_diff_bin, diff.max())
                            mean_diff_bin += diff.mean()
                            total_comparisons += 1
                    else:
                        # Fallback: compare by row index
                        min_len = min(len(df_dense), len(df_sparse_bin))
                        diff = np.abs(df_dense[col_d].iloc[:min_len] - df_sparse_bin[col_b].iloc[:min_len])
                        max_diff_bin = max(max_diff_bin, diff.max())
                        mean_diff_bin += diff.mean()
                        total_comparisons += 1
                
                if total_comparisons > 0:
                    mean_diff_bin /= total_comparisons
                
                print(f"Ground Truth (Dense) vs Sparse Binary:")
                print(f"  Max difference: {max_diff_bin:.2e}")
                print(f"  Mean difference: {mean_diff_bin:.2e}")
                if max_diff_bin < 1e-5:
                    print(f"  ✓ Results match ground truth (within numerical precision)")
                elif max_diff_bin < 1e-3:
                    print(f"  ⚠ Results differ slightly (may be due to different implementations)")
                else:
                    print(f"  ✗ Results differ significantly from ground truth")
            else:
                print(f"  ⚠ Column count mismatch: {len(score_cols_dense)} vs {len(score_cols_bin)}")
        else:
            print(f"\n✗ Sparse binary score file not found")
    except Exception as e:
        print(f"Error comparing sparse binary: {e}")
        import traceback
        traceback.print_exc()
elif not sparse_supported:
    print(f"\n⚠ Sparse binary comparison skipped (not supported)")
else:
    print(f"\n⚠ Sparse binary comparison skipped (scoring failed)")

# Compare dense vs sparse MTX (using ground truth)
if sparse_mtx_success and df_dense is not None:
    try:
        if os.path.exists(f'{bench_dir}/results_sparse_mtx.sscore'):
            # Handle header that starts with #
            with open(f'{bench_dir}/results_sparse_mtx.sscore', 'r') as f:
                header = f.readline().strip()
                if header.startswith('#'):
                    header = header[1:]
                cols_mtx = header.split('\t')
            df_sparse_mtx = pd.read_csv(f'{bench_dir}/results_sparse_mtx.sscore', sep='\t', comment='#', names=cols_mtx, skiprows=1)
            # Find score columns - match by pattern
            import re
            score_cols_mtx = [c for c in cols_mtx if '_AVG' in c and re.search(r'\d+', c)]
            def extract_num(col):
                match = re.search(r'(\d+)', col)
                return int(match.group(1)) if match else 0
            score_cols_mtx = sorted(score_cols_mtx, key=extract_num)
            
            print(f"\nSparse MTX: {len(score_cols_mtx)} score columns found")
            
            if len(score_cols_dense) == len(score_cols_mtx):
                # Match by index position (ground truth comparison)
                max_diff_mtx = 0.0
                mean_diff_mtx = 0.0
                total_comparisons = 0
                
                for idx, (col_d, col_m) in enumerate(zip(score_cols_dense, score_cols_mtx)):
                    # Match samples by IID
                    if 'IID' in df_dense.columns and 'IID' in df_sparse_mtx.columns:
                        merged = pd.merge(df_dense[['IID', col_d]], df_sparse_mtx[['IID', col_m]], on='IID', how='inner')
                        if len(merged) > 0:
                            diff = np.abs(merged[col_d] - merged[col_m])
                            max_diff_mtx = max(max_diff_mtx, diff.max())
                            mean_diff_mtx += diff.mean()
                            total_comparisons += 1
                    else:
                        # Fallback: compare by row index
                        min_len = min(len(df_dense), len(df_sparse_mtx))
                        diff = np.abs(df_dense[col_d].iloc[:min_len] - df_sparse_mtx[col_m].iloc[:min_len])
                        max_diff_mtx = max(max_diff_mtx, diff.max())
                        mean_diff_mtx += diff.mean()
                        total_comparisons += 1
                
                if total_comparisons > 0:
                    mean_diff_mtx /= total_comparisons
                
                print(f"Ground Truth (Dense) vs Sparse MTX:")
                print(f"  Max difference: {max_diff_mtx:.2e}")
                print(f"  Mean difference: {mean_diff_mtx:.2e}")
                if max_diff_mtx < 1e-5:
                    print(f"  ✓ Results match ground truth (within numerical precision)")
                elif max_diff_mtx < 1e-3:
                    print(f"  ⚠ Results differ slightly (may be due to different implementations)")
                else:
                    print(f"  ✗ Results differ significantly from ground truth")
            else:
                print(f"  ⚠ Column count mismatch: {len(score_cols_dense)} vs {len(score_cols_mtx)}")
        else:
            print(f"\n✗ Sparse MTX score file not found")
    except Exception as e:
        print(f"Error comparing sparse MTX: {e}")
        import traceback
        traceback.print_exc()
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
print(f"{'Dense (Traditional)':<25} {dense_time:<15.3f} {'1.00x':<15} {'✓':<15}")

if sparse_supported and sparse_bin_success and sparse_bin_time > 0:
    speedup = dense_time / sparse_bin_time if sparse_bin_time > 0 else 0
    speedup_str = f"{speedup:.2f}x"
    print(f"{'Sparse (Binary)':<25} {sparse_bin_time:<15.3f} {speedup_str:<15} {'✓':<15}")
elif sparse_supported:
    print(f"{'Sparse (Binary)':<25} {'N/A':<15} {'N/A':<15} {'✗ Failed':<15}")
else:
    print(f"{'Sparse (Binary)':<25} {'N/A':<15} {'N/A':<15} {'Not supported':<15}")

if sparse_supported and sparse_mtx_success and sparse_mtx_time > 0:
    speedup = dense_time / sparse_mtx_time if sparse_mtx_time > 0 else 0
    speedup_str = f"{speedup:.2f}x"
    print(f"{'Sparse (MTX)':<25} {sparse_mtx_time:<15.3f} {speedup_str:<15} {'✓':<15}")
elif sparse_supported:
    print(f"{'Sparse (MTX)':<25} {'N/A':<15} {'N/A':<15} {'✗ Failed':<15}")
else:
    print(f"{'Sparse (MTX)':<25} {'N/A':<15} {'N/A':<15} {'Not supported':<15}")

print(f"\nMemory efficiency:")
bench_dir = "${BENCH_DIR}"
actual_snps = ${ACTUAL_SNPS}
num_scores = ${NUM_SCORES}

dense_size_gb = (actual_snps * num_scores * 8) / (1024**3)
print(f"  Dense format: ~{dense_size_gb:.4f} GB")

if os.path.exists(f'{bench_dir}/sparse.binsparse'):
    sparse_size = os.path.getsize(f'{bench_dir}/sparse.binsparse')
    sparse_size_gb = sparse_size / (1024**3)
    print(f"  Sparse format: ~{sparse_size_gb:.4f} GB")
    if dense_size_gb > 0:
        print(f"  Memory savings: {(1 - sparse_size_gb/dense_size_gb)*100:.1f}%")
else:
    print(f"  Sparse format: File not found")
EOF

echo ""
echo "=========================================="
echo "  Small Test Complete"
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
    echo "Note: Sparse score matrices have been generated and are available at:"
    echo "      - ${BENCH_DIR}/sparse.mtx"
    echo "      - ${BENCH_DIR}/sparse.binsparse"
    echo "      - ${BENCH_DIR}/sparse.snp"
    echo "      - ${BENCH_DIR}/sparse.names"
fi
echo ""

