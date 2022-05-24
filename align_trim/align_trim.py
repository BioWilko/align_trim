#!/usr/bin/env python

# Written by Nick Loman (@Pathogenomenick) and packaged by Sam Wilkinson (@BioWilko)

import os
import click
from types import SimpleNamespace
import pysam

from . import scheme_downloader

from . import align_trim_funcs


@click.command()
@click.option(
    "--min-overlap",
    type=click.IntRange(min=0),
    default=0,
    help="Shortest allowed overlap to amplicon region",
)
@click.option(
    "--max-mutual-overlap",
    type=click.FloatRange(min=0, max=1),
    default=0.5,
    help="Maximum permissable proportion of overlap between read/read-pair and the next top matching amplicon as a proportion of read length",
)
@click.option("--normalise", type=int, help="Subsample to n coverage per strand")
@click.option("--report", type=click.Path(), help="Output report to file")
@click.option(
    "--trim-primers",
    default=False,
    is_flag=True,
    help="Should primers be trimmed from BAM file",
)
@click.option(
    "--verbose",
    default=False,
    is_flag=True,
    help="Print decision log messages to stderr",
)
@click.option("--prefix", type=click.STRING, default="no_prefix_supplied")
@click.option(
    "--no-read-groups",
    is_flag=True,
    default=False,
    help="Do not divide reads into groups in output",
)
@click.option("--disregard-supplemental-reads", default=True, is_flag=True)
@click.option("--debug", default=False, is_flag=True)
@click.option(
    "--enforce-amplicon-span",
    is_flag=True,
    default=False,
    help="Discard reads that do not cover the entire amplicon",
)
@click.option("--quiet", help="")
@click.option("--output_filetype", type=click.Choice(["SAM", "BAM"]), default="SAM")
@click.option(
    "--scheme-directory",
    help="Directory containing primer schemes (e.g. ~/artic/primer-schemes)",
    required=True,
)
@click.option(
    "--outdir",
    type=click.Path(exists=True),
    help="Directory in which to save outfile (supresses stdout)",
)
@click.argument("scheme")
@click.argument("infile", type=click.Path(exists=True), required=True)
def main(*_, **kwargs):
    args = SimpleNamespace(**kwargs)
    if args.infile and not os.path.exists(args.infile + ".bai"):
        pysam.index(args.infile)

    args.paired = align_trim_funcs.paired_check(args.infile)

    args.bedfile = scheme_downloader.get_scheme(args.scheme, args.scheme_directory)

    align_trim_funcs.overlap_trim(args)


if __name__ == "__main__":
    main()
