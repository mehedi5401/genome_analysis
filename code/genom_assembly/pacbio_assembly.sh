!/bin/bash -l
SBATCH -A uppmax2026-1-61
SBATCH -p pelle
SBATCH -c 4
SBATCH --mem=32G
SBATCH -t 12:00:00
SBATCH -J pacbio_assembly
SBATCH --mail-type=ALL
SBATCH --output=%x.%j.out

# Load modules
module load canu/2.3
module load SAMtools/1.22

# Define paths
WORKDIR=/home/mdha5802/1_Zhang_2017
PACBIO=$WORKDIR/genomics_data/PacBio
OUTDIR=/proj/uppmax2026-1-61/nobackup/mdha5802/assembly/pacbio

# Create output directory
mkdir -p $OUTDIR

# Step 1 — Run Canu assembly
echo "Starting Canu assembly..."
canu \
    -p E745 \
    -d $OUTDIR \
    genomeSize=2.8m \
    maxThreads=4 \
    minReadLength=1000 \
    minOverlapLength=500 \
    useGrid=false\
    -pacbio-raw $OUTDIR/*.fastq.gz

echo "Canu assembly complete"
echo "Output files:"
ls -lh $OUTDIR