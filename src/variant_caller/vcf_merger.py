#!/usr/bin/env python3
"""
Merge SNP and structural variant VCF files
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
##source=MUMmer4+CustomSVCaller
##reference={ref_name}
##contig=<ID={ref_name},length={ref_length}>
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of structural variant">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of structural variant">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
"""

def main():
    parser = argparse.ArgumentParser(description='Merge SNP and SV VCF files')
    parser.add_argument('--snps', '-s', required=True, help='Input SNPs VCF file')
    parser.add_argument('--svs', '-v', required=True, help='Input structural variants file')
    parser.add_argument('--output', '-o', required=True, help='Output merged VCF file')
    parser.add_argument('--reference', '-r', required=True, help='Reference genome file')
    
    args = parser.parse_args()
    
    try:
        ref_name, ref_length = get_reference_info(args.reference)
        logger.info(f"Reference name: {ref_name}, length: {ref_length}")
        
        with open(args.output, 'w') as out:
            out.write(write_vcf_header(ref_name, ref_length))
            
            # Read and write SNP records
            snp_records = []
            with open(args.snps, 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        # Replace old reference name with new one if needed
                        fields = line.strip().split('\t')
                        fields[0] = ref_name  # Update CHROM field
                        snp_records.append('\t'.join(fields))
            
            # Read and format SV records
            sv_records = []
            with open(args.svs, 'r') as f:
                next(f)  # Skip header
                for line in f:
                    fields = line.strip().split('\t')
                    sv_type, ref_pos, size = fields[0:3]
                    ref_pos = int(ref_pos)
                    size = int(float(size))
                    
                    if sv_type == "DEL":
                        sv_records.append(
                            f"{ref_name}\t{ref_pos}\t.\tN\t<DEL>\t.\tPASS\t"
                            f"SVTYPE=DEL;SVLEN=-{size};END={ref_pos+size}\tGT\t1/1"
                        )
                    elif sv_type == "INS":
                        sv_records.append(
                            f"{ref_name}\t{ref_pos}\t.\tN\t<INS>\t.\tPASS\t"
                            f"SVTYPE=INS;SVLEN={size}\tGT\t1/1"
                        )
                    elif sv_type == "INV":
                        sv_records.append(
                            f"{ref_name}\t{ref_pos}\t.\tN\t<INV>\t.\tPASS\t"
                            f"SVTYPE=INV;SVLEN={size};END={ref_pos+size}\tGT\t1/1"
                        )
            
            # Write sorted records
            for record in sorted(snp_records + sv_records, 
                               key=lambda x: int(x.split('\t')[1])):
                out.write(record + "\n")
                
        logger.info(f"Merged VCF file created: {args.output}")
        
    except Exception as e:
        logger.error(f"Error merging VCF files: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
    