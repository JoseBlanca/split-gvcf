import tempfile
from pathlib import Path

from test_ranges_union import VCF1, VCF2, write_compressed_file
from split_gvcf.split_genome_avoiding_vars import split_genome_avoiding_vars


def test_split_in_empty_ranges():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)
        path_vcf1 = tmp_dir_path / "1.vcf.gz"
        path_vcf2 = tmp_dir_path / "2.vcf.gz"
        write_compressed_file(path_vcf1, VCF1)
        write_compressed_file(path_vcf2, VCF2)

        chrom_size_path = tmp_dir_path / "chrom_sizes.txt"
        fhand = chrom_size_path.open("wt")
        fhand.write("20\t100\n")
        fhand.close()
        out_ranges_path = tmp_dir_path / "var_ranges.parquet"
        split_genome_avoiding_vars(
            [path_vcf1, path_vcf2],
            chrom_size_path,
            desired_segment_size=50,
            out_parquet_vars_ranges=out_ranges_path,
        )


def true_xtest():
    working_dir = Path("/vcfs/")
    vcfs_dir = working_dir
    chrom_sizes = working_dir / "chrom_sizes.txt"
    out_ranges_path = working_dir / "var_ranges.parquet"
    vcf_paths = [path for path in vcfs_dir.iterdir() if str(path).endswith(".g.vcf.gz")]
    print(vcf_paths)

    bin_working_dir = working_dir / "split_working_dir"

    split_genome_avoiding_vars(
        vcf_paths,
        bin_working_dir,
        chrom_sizes,
        desired_segment_size=1_000_000,
        out_parquet_vars_ranges=out_ranges_path,
        n_vcf_parsing_processes=6,
    )


if __name__ == "__main__":
    true_xtest()
