#!/usr/bin/env python3
"""
Convert MUMmer SNPs output to VCF format
"""
import sys
import logging
import os
from cyvcf2 import VCF
import pysam
from datetime import datetime

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_reference_info(ref_file: str):
    """Get reference sequence name and length from FASTA"""
    with pysam.FastaFile(ref_file) as fasta:
        name = fasta.references[0]  # Get first reference name
        length = fasta.get_reference_length(name)
    return name, length

def write_vcf_header(f, ref_name: str, ref_length: int):
    """Write VCF header"""
    f.write("##fileformat=VCFv4.2\n")
    f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
    f.write("##source=MUMmer4\n")
    f.write(f"##contig=<ID={ref_name},length={ref_length}>\n")
    f.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

def convert_snps_to_vcf(snps_file: str, reference_file: str, output_file: str):
    """Convert MUMmer SNPs to VCF format"""
    try:
        # Get reference info
        ref_name, ref_length = get_reference_info(reference_file)
        logger.info(f"Reference: {ref_name}, length: {ref_length}")
        
        # Process SNPs and write VCF
        with open(output_file, 'w') as vcf:
            # Write header
            write_vcf_header(vcf, ref_name, ref_length)
            
            # Process variants
            with open(snps_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                        
                    fields = line.strip().split('\t')
                    if len(fields) < 11:
                        continue
                        
                    pos = fields[0]
                    ref = fields[1]
                    alt = fields[2]
                    
                    if ref == '.' or alt == '.':
                        continue
                    
                    # Write VCF record
                    record = [
                        ref_name,          # CHROM
                        pos,               # POS
                        ".",              # ID
                        ref.upper(),      # REF
                        alt.upper(),      # ALT
                        ".",              # QUAL
                        "PASS",           # FILTER
                        "DP=30",          # INFO
                        "GT",             # FORMAT
                        "1/1"             # Sample GT
                    ]
                    vcf.write("\t".join(record) + "\n")
        
        logger.info(f"Converted SNPs to VCF: {output_file}")
        
        # Validate VCF
        try:
            vcf = VCF(output_file)
            vcf.close()
        except Exception as e:
            logger.error(f"VCF validation failed: {e}")
            raise
            
    except Exception as e:
        logger.error(f"Error converting SNPs to VCF: {e}")
        raise

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert MUMmer SNPs to VCF format")
    parser.add_argument('--input', '-i', required=True, help='Input SNPs file')
    parser.add_argument('--output', '-o', required=True, help='Output VCF file')
    parser.add_argument('--reference', '-r', required=True, help='Reference FASTA file')
    
    args = parser.parse_args()
    convert_snps_to_vcf(args.input, args.reference, args.output)

if __name__ == "__main__":
    main()
