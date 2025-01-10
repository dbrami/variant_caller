#!/usr/bin/env python3
"""
Compare variant calls between MUMmer and minimap2 aligners
"""
import argparse
import logging
from cyvcf2 import VCF
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Set, Tuple

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class VariantComparator:
    def __init__(self, mummer_vcf: str, minimap_vcf: str):
        self.mummer_vcf = mummer_vcf
        self.minimap_vcf = minimap_vcf
        
    def get_variant_set(self, vcf_file: str) -> Set[Tuple[str, int, str, str]]:
        """Get set of variants as (chrom, pos, ref, alt) tuples"""
        variants = set()
        for variant in VCF(vcf_file):
            variants.add((variant.CHROM, variant.POS, variant.REF, variant.ALT[0]))
        return variants
        
    def count_variants_by_type(self, vcf_file: str) -> Dict[str, int]:
        """Count variants by type (SNP, INS, DEL)"""
        counts = defaultdict(int)
        for variant in VCF(vcf_file):
            ref, alt = variant.REF, variant.ALT[0]
            if len(ref) == len(alt) == 1:
                counts['SNP'] += 1
            elif len(ref) > len(alt):
                counts['DEL'] += 1
            else:
                counts['INS'] += 1
        return dict(counts)
        
    def get_variant_positions(self, vcf_file: str) -> List[int]:
        """Get list of variant positions"""
        return [variant.POS for variant in VCF(vcf_file)]
        
    def plot_variant_distribution(self, mummer_pos: List[int], minimap_pos: List[int], 
                                output: str, bin_size: int = 1000):
        """Plot variant position distributions"""
        plt.figure(figsize=(12, 6))
        plt.hist([mummer_pos, minimap_pos], bins=100, label=['MUMmer', 'minimap2'], 
                alpha=0.5, density=True)
        plt.xlabel('Genome Position')
        plt.ylabel('Density')
        plt.title('Distribution of Variants Along Genome')
        plt.legend()
        plt.savefig(output)
        plt.close()
        
    def compare_variants(self) -> Dict:
        """Compare variants between the two VCF files"""
        # Get variant sets
        mummer_vars = self.get_variant_set(self.mummer_vcf)
        minimap_vars = self.get_variant_set(self.minimap_vcf)
        
        # Count variants by type
        mummer_counts = self.count_variants_by_type(self.mummer_vcf)
        minimap_counts = self.count_variants_by_type(self.minimap_vcf)
        
        # Calculate intersections
        shared_vars = mummer_vars.intersection(minimap_vars)
        mummer_only = mummer_vars - minimap_vars
        minimap_only = minimap_vars - mummer_vars
        
        # Get positions for distribution plot
        mummer_pos = self.get_variant_positions(self.mummer_vcf)
        minimap_pos = self.get_variant_positions(self.minimap_vcf)
        
        return {
            'mummer_total': len(mummer_vars),
            'minimap_total': len(minimap_vars),
            'shared': len(shared_vars),
            'mummer_only': len(mummer_only),
            'minimap_only': len(minimap_only),
            'mummer_counts': mummer_counts,
            'minimap_counts': minimap_counts,
            'positions': (mummer_pos, minimap_pos)
        }
        
    def generate_report(self, output_prefix: str):
        """Generate comparison report with plots"""
        results = self.compare_variants()
        
        # Generate report text
        with open(f"{output_prefix}_report.txt", 'w') as f:
            f.write("Variant Caller Comparison Report\n")
            f.write("==============================\n\n")
            
            f.write("Variant Counts:\n")
            f.write(f"MUMmer total variants: {results['mummer_total']}\n")
            f.write(f"minimap2 total variants: {results['minimap_total']}\n")
            f.write(f"Shared variants: {results['shared']}\n")
            f.write(f"MUMmer-only variants: {results['mummer_only']}\n")
            f.write(f"minimap2-only variants: {results['minimap_only']}\n\n")
            
            f.write("Variant Types:\n")
            f.write("MUMmer:\n")
            for vtype, count in results['mummer_counts'].items():
                f.write(f"  {vtype}: {count}\n")
            f.write("minimap2:\n")
            for vtype, count in results['minimap_counts'].items():
                f.write(f"  {vtype}: {count}\n")
        
        # Generate plots
        # Position distribution
        self.plot_variant_distribution(
            results['positions'][0], 
            results['positions'][1],
            f"{output_prefix}_distribution.png"
        )
        
        # Variant type comparison
        types = sorted(set(results['mummer_counts'].keys()) | 
                      set(results['minimap_counts'].keys()))
        comparison_data = {
            'Type': types * 2,
            'Count': [results['mummer_counts'].get(t, 0) for t in types] + 
                    [results['minimap_counts'].get(t, 0) for t in types],
            'Caller': ['MUMmer'] * len(types) + ['minimap2'] * len(types)
        }
        df = pd.DataFrame(comparison_data)
        
        plt.figure(figsize=(8, 6))
        sns.barplot(data=df, x='Type', y='Count', hue='Caller')
        plt.title('Variant Types by Caller')
        plt.savefig(f"{output_prefix}_types.png")
        plt.close()

def main():
    parser = argparse.ArgumentParser(description="Compare variant calls between aligners")
    parser.add_argument('--mummer', required=True, help='MUMmer VCF file')
    parser.add_argument('--minimap', required=True, help='minimap2 VCF file')
    parser.add_argument('--output', required=True, help='Output prefix for report and plots')
    args = parser.parse_args()
    
    comparator = VariantComparator(args.mummer, args.minimap)
    comparator.generate_report(args.output)

if __name__ == "__main__":
    main()
