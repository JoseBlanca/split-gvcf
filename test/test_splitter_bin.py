import tempfile
from pathlib import Path

from test_ranges_union import VCF1, VCF2, write_compressed_file
from split_gvcf.split_genome_avoiding_vars import split_genome_avoiding_vars


def true_xtest():
    working_dir = Path("/vcfs/")
    vcfs_dir = working_dir
    chrom_sizes = working_dir / "chrom_sizes.txt"
    vcf_paths = [path for path in vcfs_dir.iterdir() if str(path).endswith(".g.vcf.gz")]
    print(vcf_paths)

    bin_working_dir = working_dir / "split_working_dir"

    segments = split_genome_avoiding_vars(
        vcf_paths,
        bin_working_dir,
        chrom_sizes,
        desired_segment_size=1_000_000,
        n_vcf_parsing_processes=6,
    )
    for segment in segments:
        print(segment)


if __name__ == "__main__":
    true_xtest()
