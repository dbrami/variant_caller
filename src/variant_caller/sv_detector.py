#!/usr/bin/env python3
import sys
import argparse
from typing import Dict, List, Optional, Tuple
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class SVDetector:
    def __init__(self, min_size: int = 50, max_size: int = 100000):
        self.min_size = min_size
        self.max_size = max_size
        
    def parse_coords(self, coords_file: str) -> List[Dict]:
        svs = []
        try:
            with open(coords_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 13:
                        continue
                        
                    ref_start, ref_end = int(fields[0]), int(fields[1])
                    query_start, query_end = int(fields[3]), int(fields[4])
                    ref_len = abs(ref_end - ref_start)
                    query_len = abs(query_end - query_start)
                    
                    sv = self.classify_sv(ref_start, ref_end, query_start, 
                                       query_end, ref_len, query_len)
                    if sv:
                        svs.append(sv)
        except Exception as e:
            logger.error(f"Error parsing coords file: {e}")
            raise
            
        return svs
            
    def classify_sv(self, ref_start: int, ref_end: int, 
                   query_start: int, query_end: int,
                   ref_len: int, query_len: int) -> Optional[Dict]:
        size_diff = abs(ref_len - query_len)
        
        if size_diff < self.min_size or size_diff > self.max_size:
            return None
            
        if query_len > ref_len * 1.1:
            return {
                'type': 'INS',
                'ref_pos': ref_start,
                'size': query_len - ref_len,
                'ref_coords': f"{ref_start}-{ref_end}",
                'query_coords': f"{query_start}-{query_end}"
            }
        elif ref_len > query_len * 1.1:
            return {
                'type': 'DEL',
                'ref_pos': ref_start,
                'size': ref_len - query_len,
                'ref_coords': f"{ref_start}-{ref_end}",
                'query_coords': f"{query_start}-{query_end}"
            }
        elif (query_start > query_end) != (ref_start > ref_end):
            return {
                'type': 'INV',
                'ref_pos': ref_start,
                'size': ref_len,
                'ref_coords': f"{ref_start}-{ref_end}",
                'query_coords': f"{query_start}-{query_end}"
            }
        return None

def main():
    parser = argparse.ArgumentParser(description='Detect structural variants from MUMmer coords file')
    parser.add_argument('--input', '-i', required=True, help='Input coords file')
    parser.add_argument('--output', '-o', required=True, help='Output SV file')
    parser.add_argument('--min-size', type=int, default=50, help='Minimum SV size')
    parser.add_argument('--max-size', type=int, default=100000, help='Maximum SV size')
    
    args = parser.parse_args()
    
    try:
        detector = SVDetector(min_size=args.min_size, max_size=args.max_size)
        svs = detector.parse_coords(args.input)
        
        with open(args.output, 'w') as out:
            out.write("Type\tRefPosition\tSize\tRefCoords\tQueryCoords\n")
            for sv in svs:
                out.write(f"{sv['type']}\t{sv['ref_pos']}\t{sv['size']}\t"
                         f"{sv['ref_coords']}\t{sv['query_coords']}\n")
                         
        logger.info(f"Found {len(svs)} structural variants")
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
