from pathlib import Path
import functools

from genomicranges import GenomicRanges
from threaded_map_reduce import map_reduce

from split_gvcf.gvcf_parser import parse_gvcf_from_path


def unify_two_ranges(
    ranges1: GenomicRanges | None, ranges2: GenomicRanges
) -> GenomicRanges:
    if ranges1 is None:
        return ranges2
    else:
        return ranges1.union(ranges2)


def create_union_ranges_from_vcfs(
    vcf_paths: list[Path], buffer_size=10, threads=1
) -> GenomicRanges:
    if threads == 1:
        ranges = functools.reduce(
            unify_two_ranges, map(parse_gvcf_from_path, vcf_paths)
        )
    else:
        ranges = map_reduce(
            parse_gvcf_from_path,
            unify_two_ranges,
            vcf_paths,
            buffer_size=buffer_size,
            num_computing_threads=threads,
        )
    return ranges
