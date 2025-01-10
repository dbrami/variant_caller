#!/usr/bin/env python3
"""
VCF utilities using cyvcf2 for proper VCF handling
"""
import logging
from typing import Dict, List, Optional, Tuple
from cyvcf2 import VCF, Writer, Variant
import pysam
import os

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VCFHandler:
    def __init__(self, reference_file: str):
        """Initialize VCF handler with reference file"""
        self.reference_file = reference_file
        self.ref_name, self.ref_length = self._get_reference_info()
        
    def _get_reference_info(self) -> Tuple[str, int]:
        """Get reference name and length from FASTA"""
        with pysam.FastaFile(self.reference_file) as fasta:
            name = fasta.references[0]  # Get first reference name
            length = fasta.get_reference_length(name)
        return name, length

    def sort_and_normalize(self, input_vcf: str, output_vcf: str) -> None:
        """Sort and normalize VCF"""
        try:
            # Read input VCF
            vcf = VCF(input_vcf)
            
            # Collect and sort variants
            variants = []
            for variant in vcf:
                variants.append((variant.CHROM, variant.POS, variant))
            variants.sort()  # Sort by chromosome and position
            
            # Create temporary VCF with sorted variants
            temp_vcf = f"{output_vcf}.temp"
            w = Writer(temp_vcf, vcf)
            for chrom, pos, var in variants:
                w.write_record(var)
            w.close()
            vcf.close()
            
            # Replace original with sorted version
            os.rename(temp_vcf, output_vcf)
            logger.info(f"Successfully sorted and normalized VCF to {output_vcf}")
            
        except Exception as e:
            logger.error(f"Error processing VCF: {e}")
            raise

def setup_logging(level=logging.INFO):
    """Setup logging configuration"""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="VCF utilities")
    parser.add_argument('--input', '-i', help='Input VCF file')
    parser.add_argument('--output', '-o', help='Output VCF file')
    parser.add_argument('--reference', '-r', help='Reference FASTA file')
    args = parser.parse_args()
    
    handler = VCFHandler(args.reference)
    handler.sort_and_normalize(args.input, args.output)
    