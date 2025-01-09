#!/usr/bin/env python3
"""
Convert MUMmer SNPs output to VCF format
"""
import sys
import argparse
from datetime import datetime
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_reference_info(ref_file: str) -> tuple:
    """Get reference sequence name and length from FASTA"""
    ref_name = None
    length = 0
    with open(ref_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if not ref_name:  # Get name from first FASTA header
                    ref_name = line[1:].strip()  # Keep full header without '>'
            else:
                length += len(line.strip())
    return ref_name, length

def write_vcf_header(ref_name: str, ref_length: int) -> str:
    """Generate VCF header with complete contig information"""
    return f"""##fileformat=VCFv4.2
##fileDate={datetime.now().strftime('%Y%m%d')}
##source=MUMmer4
##reference={ref_name}
##contig=<ID={ref_name},length={ref_length}>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

def convert_snps_to_vcf(snps_file: str, ref_file: str, output_file: str):
    """Convert MUMmer SNPs to VCF format"""
    try:
        ref_name, ref_length = get_reference_info(ref_file)
        logger.info(f"Reference name: {ref_name}, length: {ref_length}")
        
        with open(output_file, 'w') as vcf_out:
            vcf_out.write(write_vcf_header(ref_name, ref_length))
            
            with open(snps_file, 'r') as f:
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
                        
                    vcf_record = f"{ref_name}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tDP=30\tGT\t1/1"
                    vcf_out.write(vcf_record + "\n")
                    
        logger.info(f"Converted SNPs to VCF: {output_file}")
        
    except Exception as e:
        logger.error(f"Error converting SNPs to VCF: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Convert MUMmer SNPs to VCF format')
    parser.add_argument('--input', '-i', required=True, help='Input SNPs file')
    parser.add_argument('--output', '-o', required=True, help='Output VCF file')
    parser.add_argument('--reference', '-r', required=True, help='Reference genome file')
    
    args = parser.parse_args()
    
    convert_snps_to_vcf(args.input, args.reference, args.output)

if __name__ == "__main__":
    main()
    