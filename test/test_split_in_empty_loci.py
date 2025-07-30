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
        ("chr1", 10, 25),
        ("chr1", 26, 41),
        ("chr1", 42, 49),
        ("chr1", 50, 65),
        ("chr1", 66, 81),
        ("chr1", 82, 97),
        ("chr1", 98, 100),
        ("chr2", 1, 16),
        ("chr2", 17, 32),
        ("chr2", 33, 48),
        ("chr2", 49, 59),
        ("chr2", 60, 75),
        ("chr2", 76, 91),
        ("chr2", 92, 100),
    ]
    assert list(segments) == expected
    print()
