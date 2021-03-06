#!/usr/bin/env python

from copy import copy
from collections import defaultdict, namedtuple
import sys

import numpy as np
import pysam
import pandas as pd
import click
import os

PassRead = namedtuple(
    "PassRead", ("qname", "amplicon", "coverage", "aligned_segment", "left", "right")
)
# consumesReference lookup for if a CIGAR operation consumes the reference sequence
consumesReference = [True, False, True, True, False, False, False, True]
# consumesQuery lookup for if a CIGAR operation consumes the query sequence
consumesQuery = [True, True, False, False, True, False, False, True]


def getPrimerDirection(primerID):
    """Infer the primer direction based on it's ID containing LEFT/RIGHT
    Parameters
    ----------
    primerID : string
        The primer ID from the 4th field of the primer scheme
    """
    if "LEFT" in primerID:
        return "+"
    elif "RIGHT":
        return "-"
    else:
        print("LEFT/RIGHT must be specified in Primer ID", file=sys.stderr)
        raise SystemExit(1)


def merge_sites(canonical, alt):
    """Merges a canonical primer site with an alt site, producing an interval that encompasses both
    Parameters
    ----------
    canonical : dict
        The canonical primer site, provided as a dictionary of the bed file row
    alt : dict
        The alt primer site, provided as a dictionary of the bed file row
    Returns
    -------
    dict
        A dictionary of the merged site, where the dict represents a bed file row
    """
    # base the merged site on the canonical
    mergedSite = canonical

    # check the both the canonical and alt are the same direction
    if canonical["direction"] != alt["direction"]:
        print(
            "could not merge alt with different orientation to canonical",
            file=sys.stderr,
        )
        raise SystemExit(1)

    # merge the start/ends of the alt with the canonical to get the largest window possible
    if alt["start"] < canonical["start"]:
        mergedSite["start"] = alt["start"]
    if alt["end"] > canonical["end"]:
        mergedSite["end"] = alt["end"]

    return mergedSite


def read_bed_file(fn):
    """Parses a bed file and collapses alts into canonical primer sites
    Parameters
    ----------
    fn : str
        The bedfile to parse
    Returns
    -------
    list
        A list of dictionaries, where each dictionary contains a row of the parsed bedfile.
        The available dictionary keys are - Primer_ID, direction, start, end
    """

    # read the primer scheme into a pandas dataframe and run type, length and null checks
    primers = pd.read_csv(
        fn,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "Primer_ID", "PoolName"],
        dtype={
            "chrom": str,
            "start": int,
            "end": int,
            "Primer_ID": str,
            "PoolName": str,
        },
        usecols=(0, 1, 2, 3, 4),
        skiprows=0,
    )
    if len(primers.index) < 1:
        print("primer scheme file is empty", file=sys.stderr)
        raise SystemExit(1)
    if primers.isnull().sum().sum():
        print("malformed primer scheme file", file=sys.stderr)
        raise SystemExit(1)

    # compute the direction
    primers["direction"] = primers.apply(
        lambda row: getPrimerDirection(row.Primer_ID), axis=1
    )

    # separate alt primers into a new dataframe
    altFilter = primers["Primer_ID"].str.contains("_alt")
    alts = pd.DataFrame(
        columns=("chrom", "start", "end", "Primer_ID", "PoolName", "direction")
    )
    alts = pd.concat([alts, primers[altFilter]])
    primers = primers.drop(primers[altFilter].index.values)

    # convert the primers dataframe to dictionary, indexed by Primer_ID
    #  - verify_integrity is used to prevent duplicate Primer_IDs being processed
    bedFile = primers.set_index(
        "Primer_ID", drop=False, verify_integrity=True
    ).T.to_dict()

    # if there were no alts, return the bedfile as a list of dicts
    if len(alts.index) == 0:
        return list(bedFile.values())

    # merge alts
    for _, row in alts.iterrows():
        primerID = row["Primer_ID"].split("_alt")[0]

        # check the bedFile if another version of this primer exists
        if primerID not in bedFile:

            # add to the bed file and continue
            bedFile[primerID] = row.to_dict()
            continue

        # otherwise, we've got a primer ID we've already seen so merge the alt
        mergedSite = merge_sites(bedFile[primerID], row)

        # update the bedFile
        bedFile[primerID] = mergedSite

    # return the bedFile as a list
    return [value for value in bedFile.values()]


def read_pair_generator(bam):
    """
    Generate read pairs in a BAM file.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam:
        if not read.is_proper_pair:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


def overlap_size(o1, o2):
    """
    Get the mutual overlap between two (start, end) tuples of reference co-ordinates.
    Parameters
    ----------
    o1 : tuple
        First pair of co-ordinates
    o2 : tuple
        Second pair of co-ordinates
    Returns
    -------
    int
        An integer representing the mutual overlap size between the two sets of co-ordinates.
    """
    start_1, end_1 = sorted(o1)
    start_2, end_2 = sorted(o2)
    return min(end_1, end_2) - max(start_1, start_2)


def trim(args, segment, primer_pos, end):
    """Soft mask an alignment to fit within primer start/end sites.
    Parameters
    ----------
    segment : pysam.AlignedSegment
        The aligned segment to mask
    primer_pos : int
        The position in the reference to soft mask up to (equates to the start/end position of the primer in the reference)
    end : bool
        If True, the segment is being masked from the end (i.e. for the reverse primer)
    debug : bool
        If True, will print soft masking info during trimming
    """
    # get a copy of the cigar tuples to work with
    cigar = copy(segment.cigartuples)

    # get the segment position in the reference (depends on if start or end of the segment is being processed)
    if not end:
        pos = segment.pos
    else:
        pos = segment.reference_end

    # process the CIGAR to determine how much softmasking is required
    eaten = 0
    while 1:

        # chomp CIGAR operations from the start/end of the CIGAR
        try:
            if end:
                flag, length = cigar.pop()
            else:
                flag, length = cigar.pop(0)
            if args.debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
            if args.verbose:
                print(
                    "Ran out of cigar during soft masking - completely masked read will be ignored",
                    file=sys.stderr,
                )
            break

        # if the CIGAR operation consumes the reference sequence, increment/decrement the position by the CIGAR operation length
        if consumesReference[flag]:
            if not end:
                pos += length
            else:
                pos -= length

        # if the CIGAR operation consumes the query sequence, increment the number of CIGAR operations eaten by the CIGAR operation length
        if consumesQuery[flag]:
            eaten += length

        # stop processing the CIGAR if we've gone far enough to mask the primer
        if not end and pos >= primer_pos and flag == 0:
            break
        if end and pos <= primer_pos and flag == 0:
            break

    # calculate how many extra matches are needed in the CIGAR
    extra = abs(pos - primer_pos)
    if args.debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if args.debug:
            print("Inserted a %s, %s" % (0, extra), file=sys.stderr)
        if end:
            cigar.append((0, extra))
        else:
            cigar.insert(0, (0, extra))
        eaten -= extra

    # softmask the left primer
    if not end:

        # update the position of the leftmost mappinng base
        segment.pos = pos - extra
        if args.debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if args.debug:
                print(
                    "softmask created a leading deletion in the CIGAR, shuffling the alignment",
                    file=sys.stderr,
                )
            while 1:
                if cigar[0][0] != 2:
                    break
                _, length = cigar.pop(0)
                segment.pos += length

        # now add the leading softmask
        cigar.insert(0, (4, eaten))

    # softmask the right primer
    else:
        cigar.append((4, eaten))

    # check the new CIGAR and replace the old one
    if cigar[0][1] <= 0 or cigar[-1][1] <= 0:
        raise ("invalid cigar operation created - possibly due to INDEL in primer")
    segment.cigartuples = cigar


def paired_check(bam_path):
    """
    Check the first 1000 reads of a bam file and check whether any are paired.
    """
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for i, read in enumerate(bam):
            if i >= 1000:
                break
            elif read.is_paired:
                return True
        return False


def segment_overlaps(segment, amplicons, candidate_amplicon, limits):
    overlaps = np.zeros(len(amplicons), dtype=int)

    # If candidate amplicon is not the first or last in scheme get number of overlapping positions
    # for the candidate amplicon and the amplicon either side (doing all overlaps for all reads is very slow)
    candidate_amps = (
        candidate
        for candidate in (
            candidate_amplicon - 1,
            candidate_amplicon,
            candidate_amplicon + 1,
        )
        if candidate not in limits
    )

    for amp in candidate_amps:
        candidate_record = amplicons[amp]
        overlaps[amp] = segment.get_overlap(
            candidate_record["start"], candidate_record["end"]
        )
    best, second = overlaps.argsort()[-2:][::-1]
    return overlaps, best, second


def gen_candidate_amp_map(amplicons, ref_length):
    """
    Make a mapping array of the same length as the BAM reference sequence where each position has a candidate amplicon.
    This is quite imprecise (especially near boundaries) but since we also look at the amplicon either side of the read 
    this is acceptable considering the large speedup.
    This does assume that the amplicons are named sequentially based on position, possible limitation for non-ARTIC schemes.

    TODO: Add internal amplicon numbering in cases where amplicons are not sequentially named.
    """
    map = np.zeros(ref_length + 1, dtype=int)

    for i, amp in enumerate(amplicons):
        for pos in range(amp["start"], amp["end"]):
            map[pos] = i
    return map


def secondary_segment_checks(args, amplicons, overlaps, best, second, segment):
    # don't take reads if the second best overlap is a large proportion
    # of the mutual overlap of the amplicons
    best_amplicon = amplicons[best]
    second_amplicon = amplicons[second]
    mutual = overlap_size(
        *((amp["start"], amp["end"]) for amp in (best_amplicon, second_amplicon))
    )
    if overlaps[second] > 1 / args.max_mutual_overlap * mutual:
        if args.verbose:
            print(
                "%s skipped as large secondary overlap" % (segment.qname),
                file=sys.stderr,
            )
        return False

    if (
        segment.reference_start < best_amplicon["start"]
        or segment.reference_end > best_amplicon["end"]
    ):
        if args.verbose:
            print(
                "%s skipped as read maps outside best amplicon range" % (segment.qname),
                file=sys.stderr,
            )
        return False

    return True


def first_pass(args, bam, amplicons):
    passing_reads = defaultdict(list)  # by (amplicon, is_reverse)
    n_reads = bam.get_index_statistics()[0].total / 2

    candidate_amplicon_map = gen_candidate_amp_map(amplicons, bam.lengths[0])

    limits = (
        np.amax(candidate_amplicon_map) + 1,
        np.amin(candidate_amplicon_map[np.nonzero(candidate_amplicon_map)] - 1),
    )

    if args.paired:
        # Use clicks nice progress bar so users can have a pretty progress readout
        with click.progressbar(
            read_pair_generator(bam), length=n_reads, file=sys.stderr
        ) as bam_iterator_bar:
            for segment1, segment2 in bam_iterator_bar:
                if segment_checks(args, segment1) and segment_checks(args, segment2):
                    candidate = candidate_amplicon_map[
                        (segment1.reference_start + segment2.reference_end) // 2
                    ]

                    overlaps1, __1, __2 = segment_overlaps(
                        segment1, amplicons, candidate, limits
                    )
                    overlaps2, __1, __2 = segment_overlaps(
                        segment1, amplicons, candidate, limits
                    )

                    # Combine overlap arrays within a read pair
                    overlaps = overlaps1 + overlaps2
                    best, second = overlaps.argsort()[-2:][::-1]

                else:
                    continue

                if secondary_segment_checks(
                    args, amplicons, overlaps, best, second, segment1
                ) and secondary_segment_checks(
                    args, amplicons, overlaps, best, second, segment2
                ):
                    overlap = overlaps[best]
                    amplicon = amplicons[best]

                    # check whether the insert extends to the primer at each end
                    extends = (
                        segment1.reference_start <= amplicon["insert_start"],
                        segment2.reference_end >= amplicon["insert_end"],
                    )

                    # if both primers, we call that "correctly paired"
                    if args.enforce_amplicon_span and not all(extends):
                        if args.verbose:
                            print(
                                "%s and %s skipped as pair does not fully span amplicon: %s"
                                % (segment1.query_name, segment2.query_name, extends),
                                file=sys.stderr,
                            )
                        continue
                    for segment in (segment1, segment2):
                        passing_reads[amplicon["name"], segment.is_reverse].append(
                            PassRead(
                                segment.qname,
                                amplicon["name"],
                                overlap,
                                segment,
                                *extends,
                            )
                        )

    else:
        with click.progressbar(
            bam, length=n_reads, file=sys.stderr
        ) as bam_iterator_bar:
            for segment in bam_iterator_bar:
                if segment_checks(args, segment):
                    candidate = candidate_amplicon_map[
                        (segment.reference_start + segment.reference_end) // 2
                    ]

                    overlaps, best, second = segment_overlaps(
                        segment, amplicons, candidate, limits
                    )

                else:
                    continue
                if secondary_segment_checks(
                    args, amplicons, overlaps, best, second, segment
                ):
                    overlap = overlaps[best]
                    amplicon = amplicons[best]

                    # check whether the alignment extends to the primer at each end
                    extends = (
                        segment.reference_start <= amplicon["insert_start"],
                        segment.reference_end >= amplicon["insert_end"],
                    )

                    # if both primers, we call that "correctly paired"
                    if args.enforce_amplicon_span and not all(extends):
                        if args.verbose:
                            print(
                                "%s skipped as does not span: %s"
                                % (segment.query_name, extends),
                                file=sys.stderr,
                            )
                        continue

                    if (
                        segment.reference_start < amplicon["start"]
                        or segment.reference_end > amplicon["end"]
                    ):
                        if args.verbose:
                            print(
                                "%s skipped as read-pair maps outside best amplicon range"
                                % (segment.qname),
                                file=sys.stderr,
                            )
                        continue
                    passing_reads[amplicon["name"], segment.is_reverse].append(
                        PassRead(
                            segment.qname, amplicon["name"], overlap, segment, *extends
                        )
                    )
    return passing_reads


def segment_checks(args, segment):

    if segment.is_unmapped:
        if args.verbose:
            print(
                "%s skipped as read is unmapped" % (segment.query_name),
                file=sys.stderr,
            )
        return False
    if segment.is_supplementary and args.disregard_supplemental_reads:
        if args.verbose:
            print(
                "%s skipped as supplementary" % (segment.query_name), file=sys.stderr,
            )
        return False
    return True


def generate_amplicons(args):
    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    pools = set(row["PoolName"] for row in bed)
    primer_pairs = defaultdict(dict)
    for b in bed:
        scheme, pair, side = b["Primer_ID"].split("_")
        primer_pairs[pair][side] = b
    # this data structure is more useful for searching...
    amplicons = np.fromiter(
        (
            (
                k,
                v["LEFT"]["PoolName"],
                v["LEFT"]["end"],
                v["RIGHT"]["start"],  # just insert
                v["LEFT"]["start"],
                v["RIGHT"]["end"],  # contains primers
                v["LEFT"]["Primer_ID"],
                v["RIGHT"]["Primer_ID"],
            )
            for k, v in primer_pairs.items()
        ),
        dtype=[
            ("name", int),
            ("pool", int),
            ("insert_start", int),
            ("insert_end", int),
            ("start", int),
            ("end", int),
            ("left_primer", "U20"),
            ("right_primer", "U20"),
        ],
    )
    return amplicons, pools


def overlap_trim(args):
    """Detect to which amplicon a read is derived and according to trim primers."""

    if args.report:
        reportfh = open(args.report, "w")
        reportfh.write(
            "QueryName\tReferenceStart\tReferenceEnd\t"
            "PrimerPair\t"
            "Primer1\tPrimer1Start\t"
            "Primer2\tPrimer2Start\t"
            "IsSecondary\tIsSupplementary\t"
            "Start\tEnd\tCorrectlyPaired\n"
        )

    amplicons, pools = generate_amplicons(args)

    # iterate over the alignment segments in the input BAM file
    with pysam.AlignmentFile(args.infile, "rb") as bam:
        bam_header = bam.header.copy().to_dict()
        if not args.no_read_groups:
            bam_header["RG"] = []
            for pool in pools:
                read_group = {}
                read_group["ID"] = pool
                read_group["SM"] = "pool_" + pool
                bam_header["RG"].append(read_group)

        print(
            "Reads before filtering: {}".format(bam.get_index_statistics()[0].total),
            file=sys.stderr,
        )
        passing_reads = first_pass(args, bam, amplicons)

    chosen_reads = list()

    # If
    if args.normalise:
        for (amp, is_reverse), reads in passing_reads.items():
            reads = sorted(reads, key=lambda x: x.coverage, reverse=True)
            chosen_reads.extend(reads[0 : args.normalise])
    else:
        for reads in passing_reads.values():
            chosen_reads.extend(reads)
    print("Reads after filtering: {}".format(len(chosen_reads)), file=sys.stderr)

    out_ft = "" if args.output_filetype == "SAM" else "b"

    out_extension = ".bam" if args.output_filetype == "BAM" else ".sam"

    out_path = (
        os.path.join(args.outdir, args.prefix + out_extension) if args.outdir else None
    )

    out_fh = (
        pysam.AlignmentFile("-", "w" + out_ft, header=bam_header)
        if not args.outdir
        else pysam.AlignmentFile(out_path, "w" + out_ft, header=bam_header,)
    )

    with pysam.AlignmentFile(args.infile, "rb") as bam, out_fh as outfile:
        for read in chosen_reads:
            amplicon = amplicons[amplicons["name"] == read.amplicon]

            if not args.no_read_groups:
                read.aligned_segment.set_tag("RG", str(amplicon["pool"][0]))

            if len(amplicon) > 1:
                raise IndexError(
                    "Found more than one amplicon matching: {}".format(
                        read.aligned_segment.qname
                    )
                )
            else:
                amplicon = amplicon[0]
            if not args.trim_primers:
                trim_start, trim_end = amplicon["start"], amplicon["end"]
            else:
                trim_start, trim_end = amplicon["insert_start"], amplicon["insert_end"]

            # softmask the alignment if left primer start/end inside alignment
            if read.aligned_segment.reference_start < trim_start:
                try:
                    trim(args, read.aligned_segment, trim_start, False)
                except Exception as e:
                    if args.verbose:
                        print(
                            "problem soft masking left primer in {} (error: {}), skipping".format(
                                read.aligned_segment.query_name, e
                            ),
                            file=sys.stderr,
                        )
                    continue

            # softmask the alignment if right primer start/end inside alignment
            if read.aligned_segment.reference_end > trim_end:
                try:
                    trim(args, read.aligned_segment, trim_end, True)
                except Exception as e:
                    if args.verbose:
                        print(
                            "problem soft masking right primer in {} (error: {}), skipping".format(
                                read.aligned_segment.query_name, e
                            ),
                            file=sys.stderr,
                        )
                    continue

            # check the the alignment still contains bases matching the reference
            if "M" not in read.aligned_segment.cigarstring:
                if args.verbose:
                    print(
                        "%s dropped as does not match reference post masking"
                        % (read.aligned_segment.query_name),
                        file=sys.stderr,
                    )
                continue

            # current alignment segment has passed filters, send it to the outfile
            # update the report with this alignment segment + primer details
            if args.report:
                print(
                    "QueryName\tReferenceStart\tReferenceEnd\t"
                    "PrimerPair\t"
                    "Primer1\tPrimer1Start\t"
                    "Primer2\tPrimer2Start\t"
                    "IsSecondary\tIsSupplementary\t"
                    "Start\tEnd\tCorrectlyPaired",
                    file=reportfh,
                )
                matched = (
                    1 if amplicon["left_primer"] == amplicon["right_primer"] else 0
                )
                report = "\t".join(
                    str(x)
                    for x in (
                        read.aligned_segment.query_name,
                        read.aligned_segment.reference_start,
                        read.aligned_segment.reference_end,
                        "{}_{}".format(
                            amplicon["left_primer"], amplicon["right_primer"]
                        ),
                        amplicon["left_primer"],
                        amplicon["start"],
                        amplicon["right_primer"],
                        amplicon["insert_end"],
                        read.aligned_segment.is_secondary,
                        read.aligned_segment.is_supplementary,
                        amplicon["start"],
                        amplicon["end"],
                        matched,
                    )
                )
                print(report, file=reportfh)
            outfile.write(read.segment)
    if args.outdir:
        pysam.sort(
            "-o", out_path, out_path,
        )

