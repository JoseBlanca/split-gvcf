import pytest

from genomicranges import GenomicRanges
from iranges import IRanges

from split_gvcf.ranges_operations import RangesSearcher, BeforeFirstRange, InFirstRange


def test_find_next():
    ranges = GenomicRanges(
        seqnames=["chr2", "chr2", "chr2", "chr3"],
        ranges=IRanges([100, 200, 1000, 100], [100, 100, 100, 100]),
    )

    searcher = RangesSearcher(ranges)
    with pytest.raises(BeforeFirstRange):
        searcher.find_prev_range("chr2", 90)
    with pytest.raises(InFirstRange):
        searcher.find_prev_range("chr2", 100)

    assert searcher.find_prev_range("chr2", 201) == 0
    assert searcher.find_prev_range("chr2", 301) == 1
    assert searcher.find_prev_range("chr2", 1500) == 2

    with pytest.raises(BeforeFirstRange):
        searcher.find_prev_range("chr3", 90)
    with pytest.raises(InFirstRange):
        searcher.find_prev_range("chr3", 100)

    assert searcher.find_prev_range("chr3", 1000) == 3
