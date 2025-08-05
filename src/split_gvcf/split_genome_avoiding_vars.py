from pathlib import Path
import hashlib
from multiprocessing import Pool
from subprocess import run, CalledProcessError
import os
import functools
from typing import Iterator

import polars
from genomicranges import GenomicRanges

from split_gvcf.ranges_union import unify_two_ranges
from split_gvcf.split_genome import split_in_empty_loci

VCF_PARSER_BIN = "save_var_regions_as_parquet"


def _read_chrom_sizes(chrom_size_path):
    sizes = {}
    for line in Path(chrom_size_path).open("rt"):
        chrom, size = line.split("\t")
        sizes[chrom] = int(size)
    return sizes


def _parse_gvcf(vcf_to_parse):
    cmd = [VCF_PARSER_BIN, "-i", vcf_to_parse["vcf"], "-o", vcf_to_parse["parquet"]]
    try:
        run(cmd, check=True)
    except CalledProcessError:
        if vcf_to_parse["parquet"].exists:
            os.remove(vcf_to_parse["parquet"])
        raise


def _parse_gvcfs(working_dir, vcf_paths, hashes_for_paths, n_vcf_parsing_processes):
    gvcfs_to_parse = []
    parquets = []
    for vcf in vcf_paths:
        vcf = Path(vcf)

        hash_digest = hashes_for_paths[vcf]

        parquet_path = working_dir / f"variant_regions.{hash_digest}.parquet"
        if not parquet_path.exists():
            gvcfs_to_parse.append({"vcf": vcf, "parquet": parquet_path})
        parquets.append(parquet_path)

    with Pool(processes=n_vcf_parsing_processes) as pool:
        pool.map(_parse_gvcf, gvcfs_to_parse)
    return {"parquet_paths": parquets}


def _read_vars_ranges_from_parquet(parquet: Path):
    df = polars.read_parquet(parquet)
    df = df.rename({"chrom": "seqnames", "start": "starts", "end": "ends"})
    gr = GenomicRanges.from_polars(df)
    return gr


def _create_hashes_for_paths(paths):
    hashes = {}
    for path in paths:
        path = Path(path)

        stats = path.stat()
        size_in_bytes = stats.st_size
        mtime = stats.st_mtime

        metadata_str = f"{path.name}:{size_in_bytes}:{mtime}"
        metadata_bytes = metadata_str.encode("utf-8")
        hash_digest = hashlib.sha256(metadata_bytes).hexdigest()
        hashes[path] = hash_digest
    return hashes


def _get_vars_ranges(vcf_paths, working_dir, n_vcf_parsing_processes) -> GenomicRanges:
    hashes_for_vcfs = _create_hashes_for_paths(vcf_paths)
    metadata_str = "-".join([hashes_for_vcfs[path] for path in vcf_paths])
    metadata_bytes = metadata_str.encode("utf-8")
    hash_digest = hashlib.sha256(metadata_bytes).hexdigest()
    range_union_parquet = working_dir / f"variant_regions_union.{hash_digest}.parquet"

    if not range_union_parquet.exists():
        res = _parse_gvcfs(
            working_dir, vcf_paths, hashes_for_vcfs, n_vcf_parsing_processes
        )
        ranges = functools.reduce(
            unify_two_ranges,
            map(_read_vars_ranges_from_parquet, res["parquet_paths"]),
        )
        df = polars.DataFrame(
            {
                "seqnames": ranges.get_seqnames(),
                "starts": ranges.start,
                "ends": ranges.end,
            }
        )
        df.write_parquet(range_union_parquet)
    else:
        df = polars.read_parquet(range_union_parquet)
    vars_ranges = GenomicRanges.from_polars(df)
    return vars_ranges


def create_var_ranges(vcf_paths, working_dir, n_vcf_parsing_processes) -> GenomicRanges:
    return _get_vars_ranges(vcf_paths, working_dir, n_vcf_parsing_processes)


def split_genome_avoiding_vars(
    vcf_paths: list[Path],
    working_dir: Path,
    chrom_size_path: Path,
    desired_segment_size: int,
    n_vcf_parsing_processes=1,
) -> Iterator[tuple[str, int, int]]:
    chrom_sizes = _read_chrom_sizes(chrom_size_path)

    working_dir = Path(working_dir)
    working_dir.mkdir(exist_ok=True)

    vars_ranges = _get_vars_ranges(vcf_paths, working_dir, n_vcf_parsing_processes)

    yield from split_in_empty_loci(vars_ranges, chrom_sizes, desired_segment_size)
