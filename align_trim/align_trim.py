#!/usr/bin/env python

# Written by Nick Loman (@Pathogenomenick) and packaged by Sam Wilkinson (@BioWilko)
# Thanks to Aaron Quinlan for the argparse implementation from poretools.

import argparse
import pathlib
import contextlib
import pysam
from subprocess import run
import os

from .scheme_downloader import get_scheme
from .align_trim_funcs import overlap_trim


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--temp-dir",
        help="Directory to temporarily store intermediate files (default: current dir)",
        type=pathlib.Path,
        default=os.getcwd(),
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=0,
        dest="min_overlap",
        help="Shortest allowed overlap to amplicon region",
    )
    parser.add_argument(
        "--max-mutual-overlap",
        type=int,
        default=2,
        dest="max_mutual_overlap",
        help="Maximum permissable overlap to second amplicon region as proportion of amplicon mutual overlap",
    )
    parser.add_argument(
        "--normalise", type=int, help="Subsample to n coverage per strand"
    )
    parser.add_argument("--report", type=pathlib.Path, help="Output report to file")
    parser.add_argument(
        "--trim-primers",
        action="store_true",
        help="Should primers be trimmed from BAM file",
    )
    parser.add_argument(
        "--no-read-groups",
        dest="no_read_groups",
        action="store_true",
        help="Do not divide reads into groups in output",
    )
    parser.add_argument("--debug", action="store_true", help="Debug mode")
    parser.add_argument(
        "--enforce-amplicon-span",
        action="store_true",
        dest="enforce_amplicon_span",
        help="Discard reads that do not cover the entire amplicon",
    )
    parser.add_argument(
        "--scheme-directory",
        help="Directory containing primer schemes (e.g. ~/artic/primer-schemes)",
        required=True,
    )
    parser.add_argument(
        "scheme",
        help="ARTIC primer scheme for scheme downloader (e.g. SARS-CoV-2/V4.1)",
    )
    parser.add_argument(
        "bamfile",
        help="Path to BAM file for primer-trimming / normalisation",
        type=pathlib.Path,
    )
    parser.add_argument(
        "outfile",
        help="Path to save primertrimmed / normalised BAM file",
        type=pathlib.Path,
    )
    args = parser.parse_args()
    run_pipeline(args)


def run_pipeline(args):

    if not os.path.exists(pathlib.Path.joinpath(args.bamfile, ".bai")):
        run(["samtools", "index", args.bamfile])

    args.bedfile = get_scheme(args.scheme, args.scheme_directory)

    args.tempsam = pathlib.Path.joinpath(args.temp_dir, "temp.sam")
    args.tempbam = pathlib.Path.joinpath(args.temp_dir, "temp.bam")

    overlap_trim(args)
    run(["samtools", "view", "-b", "-o", args.tempbam, args.tempsam])
    run(["samtools", "sort", "-o", args.outfile, args.tempbam])
    run(["samtools", "index", args.outfile])
    os.remove(args.tempsam)
    os.remove(args.tempbam)


if __name__ == "__main__":
    main()
