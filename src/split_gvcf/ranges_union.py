from pathlib import Path
import functools

from genomicranges import GenomicRanges


from split_gvcf.gvcf_parser import parse_gvcf_from_path


def unify_two_ranges(ranges1: GenomicRanges, ranges2: GenomicRanges) -> GenomicRanges:
    if ranges1 is None:
        return ranges2
    else:
        return ranges1.union(ranges2)


def create_union_ranges_from_vcfs(vcf_paths: list[Path]) -> GenomicRanges:
    ranges = functools.reduce(unify_two_ranges, map(parse_gvcf_from_path, vcf_paths))
    return ranges
