#!/usr/bin/env python3
"""
Convert PAF format to VCF format
"""
import sys
import argparse
from datetime import datetime
import re
import logging
from typing import List, Tuple, Dict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def get_reference_info(ref_file: str) -> tuple:
    """Get reference sequence name and length from FASTA"""
    ref_name = None
    ref_seq = ""
    with open(ref_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if not ref_name:
                    ref_name = line[1:].strip()
            else:
                ref_seq += line.strip()
    return ref_name, ref_seq

def get_reference_base(ref_seq: str, pos: int) -> str:
    """Get reference base at position (1-based position)"""
    try:
        return ref_seq[pos-1]
    except IndexError:
        logger.error(f"Position {pos} is outside reference sequence")
        return 'N'

def parse_cs_string(cs_string: str) -> List[Tuple[int, str, str]]:
    """Parse the cs string from minimap2 output"""
    ops = re.findall(r'[=\-\+\*][A-Za-z]*|\d+', cs_string)
    variants = []
    pos = 0
    
    for op in ops:
        if op[0] == '=':
            pos += len(op) - 1
        elif op[0] in '-+*':
            if op[0] == '*':
                ref = op[1]
                alt = op[2]
                variants.append((pos, ref, alt))
                pos += 1
            elif op[0] == '+':
                # For insertions, get the reference base at this position
                ref = '.'  # Will be filled in later
                alt = op[1:]
                variants.append((pos, ref, alt))
            elif op[0] == '-':
                # For deletions, get the reference sequence
                ref = op[1:]
                alt = '.'  # Will be filled in later
                variants.append((pos, ref, alt))
                pos += len(ref)
        else:
            pos += int(op)
    
    return variants

def main():
    parser = argparse.ArgumentParser(description='Convert PAF to VCF format')
    parser.add_argument('--input', '-i', required=True, help='Input PAF file')
    parser.add_argument('--output', '-o', required=True, help='Output VCF file')
    parser.add_argument('--reference', '-r', required=True, help='Reference genome file')
    
    args = parser.parse_args()
    
    try:
        # Get reference information
        ref_name, ref_seq = get_reference_info(args.reference)
        ref_length = len(ref_seq)
        logger.info(f"Reference name: {ref_name}, length: {ref_length}")
        
        with open(args.output, 'w') as vcf:
            # Write VCF header
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            vcf.write("##source=minimap2\n")
            vcf.write(f"##reference={ref_name}\n")
            vcf.write(f"##contig=<ID={ref_name},length={ref_length}>\n")
            vcf.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of variant">\n')
            vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            
            with open(args.input) as paf:
                for line in paf:
                    fields = line.strip().split("\t")
                    if len(fields) >= 12:
                        for field in fields[12:]:
                            if field.startswith("cs:Z:"):
                                cs_string = field[5:]
                                variants = parse_cs_string(cs_string)
                                
                                for pos, ref, alt in variants:
                                    actual_pos = int(fields[7]) + pos + 1
                                    
                                    # Get correct reference base(s)
                                    if ref == '.':  # Insertion
                                        ref = get_reference_base(ref_seq, actual_pos)
                                        alt = ref + alt
                                    elif alt == '.':  # Deletion
                                        ref_base = get_reference_base(ref_seq, actual_pos)
                                        alt = ref_base
                                        ref = ref_base + ref
                                    
                                    var_type = "SNP" if len(ref) == len(alt) == 1 else "INDEL"
                                    vcf.write(f"{ref_name}\t{actual_pos}\t.\t{ref}\t{alt}\t.\t"
                                            f"PASS\tTYPE={var_type}\tGT\t1/1\n")
        
        logger.info(f"Converted PAF file to VCF: {args.output}")
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

