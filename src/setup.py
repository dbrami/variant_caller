from setuptools import setup, find_packages

setup(
    name="variant_caller",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pysam",
        "biopython",
        "cyvcf2",
        "matplotlib",
        "seaborn"
    ],
    entry_points={
        'console_scripts': [
            'paf2vcf=variant_caller.paf_to_vcf:main',
            'detect_sv=variant_caller.sv_detector:main',
            'merge_vcf=variant_caller.vcf_merger:main',
            'normalize_vcf=variant_caller.vcf_utils:main',
            'compare_variants=variant_caller.comparison.compare:main'
        ],
    }
)