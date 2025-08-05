import pytest

from genomicranges import GenomicRanges
from iranges import IRanges


from split_gvcf.split_genome import split_in_empty_loci


def test_split_genome():
    ranges = GenomicRanges(
        seqnames=["chr1", "chr1", "chr2", "chr2"],
        ranges=IRanges([10, 50, 20, 60], [10, 10, 10, 10]),
    )
    chrom_sizes = {"chr1": 100, "chr2": 100}
    segments = split_in_empty_loci(
        ranges=ranges,
        chrom_sizes=chrom_sizes,
        desired_region_size=15,
    )
    expected = [
        ("chr1", 1, 9),
        ("chr1", 10, 24),
        ("chr1", 25, 39),
        ("chr1", 40, 49),
        ("chr1", 50, 64),
        ("chr1", 65, 79),
        ("chr1", 80, 94),
        ("chr1", 95, 100),
        ("chr2", 1, 15),
        ("chr2", 16, 30),
        ("chr2", 31, 45),
        ("chr2", 46, 60),
        ("chr2", 61, 75),
        ("chr2", 76, 90),
        ("chr2", 91, 100),
    ]
    assert list(segments) == expected

    chrom_sizes = {"chr1": 100, "chr2": 100, "chr3": 100}
    with pytest.raises(ValueError):
        list(
            split_in_empty_loci(
                ranges=ranges,
                chrom_sizes=chrom_sizes,
                desired_region_size=15,
            )
        )


def test_desired_larger_than_chrom():
    ranges = GenomicRanges(
        seqnames=["chr1", "chr1", "chr2", "chr2"],
        ranges=IRanges([10, 50, 20, 60], [10, 10, 10, 10]),
    )
    chrom_sizes = {"chr1": 100, "chr2": 100}
    segments = split_in_empty_loci(
        ranges=ranges,
        chrom_sizes=chrom_sizes,
        desired_region_size=200,
    )
    expected = [
        ("chr1", 1, 100),
        ("chr2", 1, 100),
    ]
    assert list(segments) == expected
