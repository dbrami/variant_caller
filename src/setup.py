from setuptools import setup, find_packages

setup(
    name="variant_caller",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas",
        "pysam",
        "biopython"
    ],
    entry_points={
        'console_scripts': [
            'paf2vcf=variant_caller.paf_to_vcf:main',
            'detect_sv=variant_caller.sv_detector:main',
            'merge_vcf=variant_caller.vcf_merger:main',
            'convert_snps=variant_caller.vcf_converter:main'
        ],
    }
)