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

## License

MIT
