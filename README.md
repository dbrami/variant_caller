# Variant Caller Pipeline

A Nextflow pipeline for variant calling using MUMmer and Minimap2.

## Prerequisites

- Nextflow
- Docker
- Python 3.9+

## Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/dbrami/variant-caller.git
   cd variant-caller
   ```

2. Build the Docker image:
   ```bash
   docker build -t macbio/variant-caller:0.1.0 -f docker/Dockerfile .
   ```

3. Run the pipeline:
   ```bash
   cd nextflow
   nextflow run main.nf --reference <reference.fa> --query <query.fa>
   ```
## Running Tests and Comparison
   ```bash
# Run both aligners
nextflow run nextflow/main.nf \
  --reference genomes/MG1655.fna \
  --query genomes/MT102.fna \
  --aligner mummer \
  --outdir test_runs/MG1655_MT102_mummer

nextflow run nextflow/main.nf \
  --reference genomes/MG1655.fna \
  --query genomes/MT102.fna \
  --aligner minimap2 \
  --outdir test_runs/MG1655_MT102_minimap2

# compare results
  compare_variants \
  --mummer test_runs/MG1655_MT102_mummer/final/normalized.vcf.gz \
  --minimap test_runs/MG1655_MT102_minimap2/final/normalized.vcf.gz \
  --output test_runs/comparison
   ```

## License

MIT
