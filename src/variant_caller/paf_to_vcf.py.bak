#!/usr/bin/env python3
"""
Convert PAF format to VCF format
"""
import sys
import argparse
from datetime import datetime
import re
import logging
from typing import List, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
                ref = 'N'
                alt = op[1:]
                variants.append((pos, ref, alt))
            elif op[0] == '-':
                ref = op[1:]
                alt = 'N'
                variants.append((pos, ref, alt))
                pos += len(ref)
        else:
            pos += int(op)
    
    return variants

def main():
    parser = argparse.ArgumentParser(description='Convert PAF to VCF format')
    parser.add_argument('--input', '-i', required=True, help='Input PAF file')
    parser.add_argument('--output', '-o', required=True, help='Output VCF file')
    parser.add_argument('--reference', '-r', required=True, help='Reference genome name')
    
    args = parser.parse_args()
    
    try:
        with open(args.output, 'w') as vcf:
            # Write VCF header
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            vcf.write("##source=minimap2\n")
            vcf.write(f"##reference={args.reference}\n")
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
                                    var_type = "SNP" if len(ref) == len(alt) == 1 else "INDEL"
                                    vcf.write(f"{args.reference}\t{actual_pos}\t.\t{ref}\t{alt}\t.\t"
                                            f"PASS\tTYPE={var_type}\tGT\t1/1\n")
        
        logger.info(f"Converted PAF file to VCF: {args.output}")
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
