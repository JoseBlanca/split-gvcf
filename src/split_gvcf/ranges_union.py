from pathlib import Path
import functools

from genomicranges import GenomicRanges


def unify_two_ranges(ranges1: GenomicRanges, ranges2: GenomicRanges) -> GenomicRanges:
    if ranges1 is None:
        return ranges2
    else:
        return ranges1.union(ranges2)
