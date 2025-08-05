import tempfile
from pathlib import Path

from test_ranges_union import write_compressed_file
from split_gvcf.split_genome_avoiding_vars import split_genome_avoiding_vars

VCF1 = b"""##fileformat=VCFv4.5
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t10\trs6054257\tAG\tA,<NON_REF>\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t20\t.\tT\tT\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t30\t.\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t40\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t50\tindel\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
21\t10\tindel\tGTC\tG,GTC\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
"""

VCF2 = b"""##fileformat=VCFv4.5
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
20\t10\trs6054257\tGAA\tA,<NON_REF>\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|2:48:1:51,51\t3|4:48:8:51,51\t5/6000:43:5:.,.
20\t20\t.\tT\t<NON_REF>\t.\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t.|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
20\t29\t.\tTA\tT,TAC,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
20\t40\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
20\t50\tindel\tGTC\tG,GTCT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t0/1:35:4\t0/2:17:2\t1/1:40:3
"""


def test_split_genome():
    with tempfile.TemporaryDirectory() as working_dir:
        working_dir_path = Path(working_dir)
        vcf1 = working_dir_path / "1.g.vcf.gz"
        write_compressed_file(vcf1, VCF1)
        vcf2 = working_dir_path / "2.g.vcf.gz"
        write_compressed_file(vcf2, VCF2)
        chrom_sizes_path = working_dir_path / "chrom_sizes.txt"
        chrom_sizes_fhand = chrom_sizes_path.open("wt")
        chrom_sizes_fhand.write("20\t65\n21\t20\n")
        chrom_sizes_fhand.flush()

        segments = list(
            split_genome_avoiding_vars(
                vcf_paths=[vcf1, vcf2],
                working_dir=working_dir_path,
                chrom_size_path=chrom_sizes_path,
                desired_segment_size=9,
            )
        )
        expected = [
            ("20", 1, 9),
            ("20", 10, 18),
            ("20", 19, 27),
            ("20", 28, 36),
            ("20", 37, 45),
            ("20", 46, 54),
            ("20", 55, 63),
            ("20", 64, 65),
            ("21", 1, 9),
            ("21", 10, 18),
            ("21", 19, 20),
        ]
        assert segments == expected
        segments = list(
            split_genome_avoiding_vars(
                vcf_paths=[vcf1, vcf2],
                working_dir=working_dir_path,
                chrom_size_path=chrom_sizes_path,
                desired_segment_size=10,
            )
        )
        expected = [
            ("20", 1, 9),
            ("20", 10, 19),
            ("20", 20, 29),
            ("20", 30, 39),
            ("20", 40, 49),
            ("20", 50, 59),
            ("20", 60, 65),
            ("21", 1, 9),
            ("21", 10, 19),
            ("21", 20, 20),
        ]
        assert segments == expected
        segments = list(
            split_genome_avoiding_vars(
                vcf_paths=[vcf1, vcf2],
                working_dir=working_dir_path,
                chrom_size_path=chrom_sizes_path,
                desired_segment_size=11,
            )
        )
        expected = [
            ("20", 1, 9),
            ("20", 10, 20),
            ("20", 21, 31),
            ("20", 32, 42),
            ("20", 43, 53),
            ("20", 54, 64),
            ("20", 65, 65),
            ("21", 1, 9),
            ("21", 10, 20),
            ("21", 21, 20),
        ]
        assert segments == expected
        segments = list(
            split_genome_avoiding_vars(
                vcf_paths=[vcf1, vcf2],
                working_dir=working_dir_path,
                chrom_size_path=chrom_sizes_path,
                desired_segment_size=12,
            )
        )
        expected = [
            ("20", 1, 9),
            ("20", 10, 21),
            ("20", 22, 33),
            ("20", 34, 45),
            ("20", 46, 57),
            ("20", 58, 65),
            ("21", 1, 9),
            ("21", 10, 20),
        ]
        assert segments == expected
        segments = list(
            split_genome_avoiding_vars(
                vcf_paths=[vcf1, vcf2],
                working_dir=working_dir_path,
                chrom_size_path=chrom_sizes_path,
                desired_segment_size=13,
            )
        )
        expected = [
            ("20", 1, 13),
            ("20", 14, 26),
            ("20", 27, 39),
            ("20", 40, 49),
            ("20", 50, 62),
            ("20", 63, 65),
            ("21", 1, 13),
            ("21", 14, 20),
        ]
        assert segments == expected
