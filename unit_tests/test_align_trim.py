from align_trim import align_trim_funcs
import os
import pysam
from types import SimpleNamespace
import numpy as np
import pandas as pd

from pandas.testing import assert_frame_equal

THIS_DIR = os.path.dirname(os.path.realpath(__file__))

args = SimpleNamespace(
    bedfile=os.path.join(THIS_DIR, "test_bed_file.bed"),
    infile=os.path.join(THIS_DIR, "test_bam_single.bam"),
    paired=False,
    max_mutual_overlap=0.5,
    min_overlap=0,
    trim_primers=True,
    discard_supplemental_reads=True,
    enforce_amplicon_span=True,
    debug=False,
)


class TestClass:
    def test_getPrimerDirection(self):
        primerID = "nCoV-2019_5_LEFT"
        assert align_trim_funcs.getPrimerDirection(primerID) == "+"

    def test_merge_sites(self):
        canonical = {
            "chrom": "MN908947.3",
            "Primer_ID": "nCoV-2019_21_LEFT",
            "direction": "+",
            "start": "6167",
            "end": "6196",
        }
        alt = {
            "chrom": "MN908947.3",
            "Primer_ID": "nCoV-2019_21_LEFT_alt2",
            "direction": "+",
            "start": "6168",
            "end": "6197",
        }
        expected_output = {
            "chrom": "MN908947.3",
            "Primer_ID": "nCoV-2019_21_LEFT",
            "direction": "+",
            "start": "6167",
            "end": "6197",
        }
        mergedSite = align_trim_funcs.merge_sites(canonical, alt)
        assert mergedSite == expected_output

    def test_read_bed_file(self):
        test_path = os.path.join(THIS_DIR, "test_bed_file.bed")
        expected_bed_path = os.path.join(THIS_DIR, "expected_read_bed_file_output.csv",)

        expected_output = pd.read_csv(
            expected_bed_path, keep_default_na=False, dtype={"PoolName": "object"}
        )

        test_output = align_trim_funcs.read_bed_file(test_path)
        test_output = pd.DataFrame(test_output)

        assert_frame_equal(expected_output, test_output)

        # check_dtype=False,
        # check_frame_type=False,
        # check_column_type=False,
        # check_names=False,

    def test_overlap_size(self):
        o1 = (500, 520)
        o2 = (550, 510)

        output = align_trim_funcs.overlap_size(o1, o2)
        assert output == 10

    def test_trim(self):
        expected_out = "9366049c-71fe-4cc5-a31f-613ba156133b	16	MN908947.3	28850	60	74S18M1D64M1I32M1I57M1D1M1D1M1D5M1D14M1D26M2I3M1I5M2D1M1D10M1D31M1D39M1D59M69S	*	0	0	TTACGGTATTGCTAAGGTTAATAGGGAAACACGATAGAATCCGAACAGCACCTTCCTCATCACGTAGTCGCAACAGTTCAAGAAATTCAACTCAGGCAGCAGTATGGGAACTTCTCCTGCTAGAATGGCTGGCAATGGCTGTGATGCTGCTCTTGCTTTTGCTGCTGCTTGACAGATTGAACCAGCTTGCAGAGCAAAATGTCTGGTAAAGGCCAACAACAACAAGGCCAAACTGTCACTAAGAAATTCGCCGAGTTTCTAAGAAGCCCGGCAAAAACGTACTGCCACTAAAGCTAATAGCAATGAATATGCTTTCGCAGACGTGGTCCAGAACAAACCCAAGGAAATTTGGGGACCAGGAACTAATCAGACAAGGAACTGATTACAACATTGGCCGCAAATTGCACAATTTGCCCCCAGCGCTTCAGCGTTCTTCGGAATGTCGAGGTGCTGTTCGGATTCTATCGTGTTTCCCTCATTAACCTTAGCAATACGTAACTGAACGAAGCACATT	(*+,1-,-64=;3@A;,;:32,0=@@DFFBF?>E@<9>=;:8?AC??@09/49;7:9/<**-:79=C?:@FHE@?=DFAA=<??='*,:88++*/87;96--&&&073>?8@84?<B./17:;>AC>='@-?5<;;==,>0;:<;2,5/8;9CC;&'7EE6487296;84-+9=:=;98=;?==C:?>;5>;?9=@CB98:=618?AB<<<9?@A>88707<71;8787;D-7*.87'54985;<;8('1%$&()$*+*(6:7855-'%45:CFKJIEC>@<8:??<98<75$$%&$%($$('&*1.,)((,..1-.*=<@=DBJKA>BDBFCA?L@B6>CC::CA@AII@??CD3.9:+?B><DDDDI>CAGC;<CA?;8),886@E?<<:DD=;B8<939:87E7+86;:)9>:::;CD=>B>748>A:;;=G:+*EC?>=8:;::111*.+.-;;A?><76.-3975AC766-+7566756,.832AA?C>CA<+;/3D<:5*+*')%%%%	NM:i:23	ms:i:646	AS:i:646	nn:i:0	tp:A:P	cm:i:37	s1:i:266	s2:i:0	de:f:0.0521	rl:i:0"

        bam_path = os.path.join(THIS_DIR, "test_bam_single.bam")
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for segment in bam:
                align_trim_funcs.trim(args, segment, 28827, False)
                align_trim_funcs.trim(args, segment, 28849, False)
                post_trimming_segment = segment.to_string()

        assert post_trimming_segment == expected_out

    def test_paired_check(self):
        single_path = os.path.join(THIS_DIR, "test_bam_single.bam")
        paired_path = os.path.join(THIS_DIR, "test_bam_paired.bam")
        assert align_trim_funcs.paired_check(single_path) == False
        assert align_trim_funcs.paired_check(paired_path) == True

    def test_segment_overlaps(self):
        args = SimpleNamespace(bedfile=os.path.join(THIS_DIR, "test_bed_file.bed"))
        single_path = os.path.join(THIS_DIR, "test_bam_single.bam")

        with pysam.AlignmentFile(single_path, "rb") as bam:
            for segment in bam:
                if segment.query_name == "9366049c-71fe-4cc5-a31f-613ba156133b":
                    break

        amplicons, _ = align_trim_funcs.generate_amplicons(args)

        expected_output = np.zeros(len(amplicons), dtype=int)
        expected_output[95] = 86
        expected_output[96] = 388
        expected_output[97] = 90
        expected_output = expected_output.tolist()

        out_overlaps, out_best, out_second = align_trim_funcs.segment_overlaps(
            segment, amplicons, 96, (0, 99)
        )

        out_overlaps = out_overlaps.tolist()

        assert out_overlaps == expected_output
        assert out_best == 96
        assert out_second == 97

    def test_generate_amplicons(self):
        args = SimpleNamespace(bedfile=os.path.join(THIS_DIR, "test_bed_file.bed"))

        amplicon_csv_path = os.path.join(THIS_DIR, "expected_amplicons.csv")

        with open(amplicon_csv_path, "rt") as amplicon_csv_fh:
            expected_amplicons = pd.read_csv(amplicon_csv_fh, keep_default_na=False)

        amplicons, pools = align_trim_funcs.generate_amplicons(args)

        amplicons = pd.DataFrame(amplicons)

        assert_frame_equal(expected_amplicons, amplicons)

        assert pools == {"1", "2"}

    def test_first_pass(self):
        ### Add reads to the test bam which should be thrown out for each possible reason, check that they are with verbose
        amplicons, _ = align_trim_funcs.generate_amplicons(args)

        # expected_output =
        with pysam.AlignmentFile(args.infile, "rb") as bam:

            passing_reads = align_trim_funcs.first_pass(args, bam, amplicons)
            print(passing_reads)
        assert passing_reads == expected_output
