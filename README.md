# Variant Calling Pipeline

A Nextflow pipeline for variant calling between two complete genomes using MUMmer or Minimap2.

## Features

- SNP and small indel detection
- Optional structural variant (SV) detection
- Support for both MUMmer and Minimap2 aligners
- Standardized VCF output
- Automatic sorting and normalization of variants

## Prerequisites

- Nextflow
- Docker
- Python 3.9+

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/variant-caller.git
cd variant-caller

# Build Docker image
docker build -t macbio/variant-caller:0.1.0 -f docker/Dockerfile .
```

## Usage

Basic SNP calling:
```bash
nextflow run main.nf \
  --reference ref.fa \
  --query query.fa
```

Include structural variant detection:
```bash
nextflow run main.nf \
  --reference ref.fa \
  --query query.fa \
  --detect_sv
```

Use Minimap2 instead of MUMmer:
```bash
nextflow run main.nf \
  --reference ref.fa \
  --query query.fa \
  --aligner minimap2
```

## Parameters

### Required Arguments
- `--reference`: Path to reference genome (FASTA format)
- `--query`: Path to query genome (FASTA format)

### Optional Arguments
- `--aligner`: Choice of alignment tool (mummer/minimap2) [default: mummer]
- `--outdir`: Output directory [default: results]
- `--detect_sv`: Enable structural variant detection [default: false]

### SV Detection Options
These parameters are only used when `--detect_sv` is enabled:
- `--min_sv_size`: Minimum size for structural variants [default: 50]
- `--max_sv_size`: Maximum size for structural variants [default: 100000]

## Output

The pipeline generates:
- Alignment files (delta/PAF format)
- Variant calls in VCF format
- Normalized and sorted final VCF
- Optional structural variant calls

## Docker Containers

The pipeline uses the following Docker containers:
- `macbio/mummer-arm64`: MUMmer aligner
- `macbio/minimap2-arm64`: Minimap2 aligner
- `macbio/variant-caller`: Custom variant processing tools

## License

MIT
