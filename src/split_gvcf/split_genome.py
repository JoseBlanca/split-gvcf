from genomicranges import GenomicRanges

from split_gvcf.ranges_operations import RangesSearcher, BeforeFirstRange, InFirstRange


def split_in_empty_loci(
    ranges: GenomicRanges,
    chrom_sizes: dict[str, int],
    desired_region_size: int,
):
    searcher = RangesSearcher(ranges)

    if set(chrom_sizes.keys()).difference(searcher.sorted_seqnames):
        raise ValueError(
            "There are chromosome with sizes that are not found in the ranges"
        )

    for chrom in searcher.sorted_seqnames:
        segment_start = 1

        while True:
            segment_end = segment_start + desired_region_size

            prev_segment = None
            try:
                prev_segment = searcher.find_prev_range(chrom, segment_end)
            except BeforeFirstRange:
                pass
            except InFirstRange:
                idx = searcher.first_range_for_seq[chrom]
                segment_end = int(searcher.starts[idx]) - 1

            if prev_segment is not None:
                if prev_segment < len(ranges) - 1:
                    next_segment = prev_segment + 1
                    if (
                        searcher.starts[next_segment] < segment_end
                        and chrom == searcher.seqnames[next_segment]
                    ):
                        segment_end = int(searcher.starts[next_segment]) - 1

            chrom_done = False
            if segment_end > chrom_sizes[chrom]:
                segment_end = chrom_sizes[chrom]
                chrom_done = True

            yield chrom, segment_start, segment_end
            if chrom_done:
                segment_start = 1
                break

            segment_start = segment_end + 1
