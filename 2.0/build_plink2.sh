#!/bin/bash
source /n/home03/btruong/miniconda3_231202/etc/profile.d/conda.sh
conda activate polyads
cd /n/holystore01/LABS/price_lab/Lab/btruong/Tools/plink-largescore1/2.0/

echo "Cleaning..."
make clean > build_script.log 2>&1

echo "Building..."
# Serial build to avoid race conditions and ensure clear logging
make BASEFLAGS='-DZSTD_MULTITHREAD -DIGNORE_BUNDLED_ZSTD -DNOLAPACK' BLASFLAGS64='' >> build_script.log 2>&1

echo "Build finished. Checking binary..."
ls -l bin/plink2 >> build_script.log 2>&1
if [ -f bin/plink2 ]; then
    echo "SUCCESS: bin/plink2 exists." >> build_script.log
else
    echo "FAILURE: bin/plink2 does not exist." >> build_script.log
fi
