from enum import Enum
from typing import Iterator, Tuple
from array import array
import gzip

import iranges
from genomicranges import GenomicRanges

FILTER_FIELD_IDX = 5
NON_REF = b"<NON_REF>"


def file_is_gzipped(fhand):
    magic = fhand.read(2)
    fhand.seek(0)
    return magic == b"\x1f\x8b"


class VcfSection(Enum):
    HEADER = 1
    BODY = 2


def parse_gvcf(fhand) -> Iterator[Tuple[str, int, int]]:
    if file_is_gzipped(fhand):
        fhand = gzip.GzipFile(fileobj=fhand)

    section = VcfSection.HEADER
    for line in fhand:
        if section == VcfSection.HEADER:
            if line.startswith(b"##"):
                continue
            elif line.startswith(b"#CHROM"):
                section = VcfSection.BODY
                continue
            else:
                raise RuntimeError(
                    "In section body lines that not start with # are not allowed"
                )

        chrom, pos, _, ref, alt, _ = line.split(b"\t", FILTER_FIELD_IDX)
        chrom = chrom.decode()
        pos = int(pos)
        if alt == NON_REF:
            continue
        alt_alleles = filter(lambda x: x != NON_REF, alt.split(b","))
        max_alt_allele_len = max(map(len, alt_alleles))
        max_allele_len = max((len(ref), max_alt_allele_len))
        yield chrom, pos, max_allele_len


def parse_gvcf_into_ranges(fhand) -> GenomicRanges:
    seq_names = []
    starts = array("L")
    widths = array("L")
    for var in parse_gvcf(fhand):
        seq_names.append(var[0])
        _, pos, max_allele_len = var
        starts.append(pos)
        widths.append(max_allele_len)
    ranges = iranges.IRanges(starts, widths)
    ranges = GenomicRanges(seqnames=seq_names, ranges=ranges)
    return ranges
