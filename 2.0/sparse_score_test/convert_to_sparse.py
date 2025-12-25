
import sys
import os

def convert_to_sparse(dense_file, out_prefix):
    # Read dense file
    # Format: ID Allele Effect1 Effect2 ...
    # Skip header
    # Sparse format:
    # .mtx: Matrix Market
    # .snp: ID Allele
    # .names: Score names

    with open(dense_file, 'r') as f:
        header_line = f.readline().strip()
        header = header_line.split()
        score_names = header[2:]
        
        snps = []
        entries = [] # (row_idx, col_idx, val)
        
        row_idx = 1
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3: continue
            
            sid = parts[0]
            allele = parts[1]
            snps.append((sid, allele))
            
            for col_idx, val_str in enumerate(parts[2:]):
                try:
                   val = float(val_str)
                   if val != 0.0:
                       entries.append((row_idx, col_idx + 1, val))
                except ValueError:
                   continue
            row_idx += 1
            
    num_rows = len(snps)
    num_cols = len(score_names)
    num_entries = len(entries)
    
    # Write .mtx
    with open(f"{out_prefix}.mtx", 'w') as f:
        f.write("%%MatrixMarket matrix coordinate real general\n")
        f.write(f"{num_rows} {num_cols} {num_entries}\n")
        for r, c, v in entries:
            f.write(f"{r} {c} {v}\n")
            
    # Write .snp
    with open(f"{out_prefix}.snp", 'w') as f:
        for sid, allele in snps:
            f.write(f"{sid}\t{allele}\n")
            
    # Write .names
    with open(f"{out_prefix}.names", 'w') as f:
        for name in score_names:
            f.write(f"{name}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_sparse.py <dense_file> <out_prefix>")
        sys.exit(1)
    
    convert_to_sparse(sys.argv[1], sys.argv[2])
