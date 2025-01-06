#!/bin/bash

# Verification script for variant-caller setup
set -e

echo "Verifying variant-caller setup..."

# Check directory structure
echo "Checking directory structure..."
for dir in nextflow src/variant_caller src/tests docker; do
    if [ ! -d "$dir" ]; then
        echo "ERROR: Directory $dir not found!"
        exit 1
    fi
done
echo "✓ Directory structure OK"

# Check required files exist
echo "Checking required files..."
required_files=(
    "src/setup.py"
    "src/requirements.txt"
    "src/variant_caller/sv_detector.py"
    "src/variant_caller/paf_to_vcf.py"
    "src/variant_caller/vcf_merger.py"
    "src/variant_caller/__init__.py"
    "docker/Dockerfile"
    "nextflow/nextflow.config"
    "nextflow/containers.config"
    "README.md"
)

for file in "${required_files[@]}"; do
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file $file not found!"
        exit 1
    fi
done
echo "✓ Required files OK"

# Check Python package installation
echo "Checking Python package..."
if ! python3 -m pip install -e src/; then
    echo "ERROR: Failed to install Python package!"
    exit 1
fi
echo "✓ Python package installation OK"

# Check if commands are available
echo "Checking command availability..."
commands=("paf2vcf" "detect_sv" "merge_vcf")
for cmd in "${commands[@]}"; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Command $cmd not found in PATH!"
        exit 1
    fi
done
echo "✓ Commands available OK"

# Check Docker
echo "Checking Docker..."
if ! docker info >/dev/null 2>&1; then
    echo "ERROR: Docker is not running or not accessible!"
    exit 1
fi
echo "✓ Docker OK"

# Check Docker image
echo "Checking Docker image..."
if ! docker images | grep -q "variant-caller"; then
    echo "ERROR: Docker image variant-caller not found!"
    exit 1
fi
echo "✓ Docker image OK"

echo "All checks passed successfully!"
echo 
echo "Next steps:"
echo "1. Push Docker image: docker push macbio/variant-caller:0.1.0"
echo "2. Create GitHub repository and push code"
echo "3. Test pipeline with sample data"
