#!/usr/bin/env python

from copy import copy
from collections import defaultdict, namedtuple
import sys

import numpy as np
import pysam
import pandas as pd


PassRead = namedtuple("PassRead", ("qname", "amplicon", "coverage", "left", "right"))
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


def overlap_size(o1, o2):
    s1, e1 = sorted(o1)
    s2, e2 = sorted(o2)
    return min(e1, e2) - max(s1, s2)


def trim(segment, primer_pos, end, debug):
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
            if debug:
                print("Chomped a %s, %s" % (flag, length), file=sys.stderr)
        except IndexError:
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
    if debug:
        print("extra %s" % (extra), file=sys.stderr)
    if extra:
        if debug:
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
        if debug:
            print("New pos: %s" % (segment.pos), file=sys.stderr)

        # if proposed softmask leads straight into a deletion, shuffle leftmost mapping base along and ignore the deletion
        if cigar[0][0] == 2:
            if debug:
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
    return


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
            "Start\tEnd\tCorrectlyPaired"
        )

    # open the primer scheme and get the pools
    bed = read_bed_file(args.bedfile)
    pools = set(row["PoolName"] for row in bed)
    pools.add("unmatched")
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
            ("pool", str),
            ("insert_start", int),
            ("insert_end", int),
            ("start", int),
            ("end", int),
            ("left_primer", "U20"),
            ("right_primer", "U20"),
        ],
    )

    # iterate over the alignment segments in the input SAM file
    passing_reads = defaultdict(list)  # by (amplicon, is_reverse)
    with pysam.AlignmentFile(args.bamfile, "rb") as bam:
        bam_header = bam.header.copy().to_dict()
        if not args.no_read_groups:
            bam_header["RG"] = []
            for pool in pools:
                read_group = {}
                read_group["ID"] = pool
                bam_header["RG"].append(read_group)

        for segment in bam:
            # filter out unmapped and supplementary alignment segments
            if segment.is_unmapped:
                print("%s skipped as unmapped" % (segment.query_name), file=sys.stderr)
                continue
            if segment.is_supplementary:
                print(
                    "%s skipped as supplementary" % (segment.query_name),
                    file=sys.stderr,
                )
                continue

            # determine the amplicon by largest overlap (and find second best)
            overlaps = np.zeros(len(amplicons), dtype=int)
            for i, amp in enumerate(amplicons):
                overlaps[i] = segment.get_overlap(amp["start"], amp["end"])
            best, second = np.argpartition(overlaps, -2)[:-3:-1]
            if overlaps[best] < args.min_overlap:
                print("%s skipped as no good overlap" % segment.qname, file=sys.stderr)
                continue

            # don't take reads if the second best overlap is a large proportion
            # of the mutual overlap of the amplicons
            best_amplicon = amplicons[best]
            second_amplicon = amplicons[second]
            mutual = overlap_size(
                *(
                    (amp["start"], amp["end"])
                    for amp in (best_amplicon, second_amplicon)
                )
            )
            if overlaps[second] > args.max_mutual_overlap * mutual:
                print(
                    "%s skipped as large secondary overlap" % segment.qname,
                    file=sys.stderr,
                )
                continue
            overlap = overlaps[best]
            amplicon = amplicons[best]

            # check whether the alignment extends to the primer at each end
            extends = (
                segment.reference_start <= amplicon["insert_start"],
                segment.reference_end >= amplicon["insert_end"],
            )

            # if both primers, we call that "correctly paired"
            if args.enforce_amplicon_span and not all(extends):
                print(
                    "%s skipped as does not span: %s" % (segment.query_name, extends),
                    file=sys.stderr,
                )
                continue
            passing_reads[amplicon["name"], segment.is_reverse].append(
                PassRead(segment.qname, amplicon["name"], overlap, *extends)
            )
    # end first pass

    # filter alignments
    print(
        "Reads before filtering: {}".format(
            sum(len(x) for x in passing_reads.values())
        ),
        file=sys.stderr,
    )
    chosen_reads = list()
    if args.normalise:
        for (amp, is_reverse), reads in passing_reads.items():
            reads = sorted(reads, key=lambda x: x.coverage, reverse=True)
            chosen_reads.extend(reads[0 : args.normalise])
    else:
        for reads in passing_reads.values():
            chosen_reads.extend(reads)
    print("Reads after filtering: {}".format(len(chosen_reads)), file=sys.stderr)
    chosen_reads = {r.qname: r for r in chosen_reads}

    with pysam.AlignmentFile(args.bamfile, "rb") as bam, pysam.AlignmentFile(
        args.tempsam, "wh", header=bam_header
    ) as outfile:
        for segment in bam:
            wanted = segment.qname in chosen_reads
            if not wanted or segment.is_unmapped or segment.is_supplementary:
                continue
            chosen = chosen_reads[segment.qname]
            amplicon = amplicons[amplicons["name"] == chosen.amplicon]
            if not args.no_read_groups:
                segment.set_tag("RG", str(amplicon["pool"][0]))

            if len(amplicon) > 1:
                raise IndexError(
                    "Found more than one amplicon matching: {}".format(chosen)
                )
            else:
                amplicon = amplicon[0]
            if args.trim_primers:
                trim_start, trim_end = amplicon["start"], amplicon["end"]
            else:
                trim_start, trim_end = amplicon["insert_start"], amplicon["insert_end"]

            # softmask the alignment if left primer start/end inside alignment
            if segment.reference_start < trim_start:
                try:
                    trim(segment, trim_start, False, args.debug)
                except Exception as e:
                    print(
                        "problem soft masking left primer in {} (error: {}), skipping".format(
                            segment.query_name, e
                        ),
                        file=sys.stderr,
                    )
                    continue

            # softmask the alignment if right primer start/end inside alignment
            if segment.reference_end > trim_end:
                try:
                    trim(segment, trim_end, True, args.debug)
                except Exception as e:
                    print(
                        "problem soft masking right primer in {} (error: {}), skipping".format(
                            segment.query_name, e
                        ),
                        file=sys.stderr,
                    )
                    continue

            # check the the alignment still contains bases matching the reference
            if "M" not in segment.cigarstring:
                print(
                    "%s dropped as does not match reference post masking"
                    % (segment.query_name),
                    file=sys.stderr,
                )
                continue

            # current alignment segment has passed filters, send it to the outfile
            # update the report with this alignment segment + primer details
            if args.report:
                # "QueryName\tReferenceStart\tReferenceEnd\t"
                # "PrimerPair\t"
                # "Primer1\tPrimer1Start\t"
                # "Primer2\tPrimer2Start\t"
                # "IsSecondary\tIsSupplementary\t"
                # "Start\tEnd\tCorrectlyPaired", file=reportfh)
                report = "\t".join(
                    str(x)
                    for x in (
                        segment.query_name,
                        segment.reference_start,
                        segment.reference_end,
                        "{}_{}".format(
                            amplicon["left_primer"], amplicon["right_primer"]
                        ),
                        amplicon["left_primer"],
                        amplicon["start"],
                        amplicon["right_primer"],
                        amplicon["insert_end"],
                        segment.is_secondary,
                        segment.is_supplementary,
                        amplicon["start"],
                        amplicon["end"],
                        int(all(extends)),
                    )
                )
                if args.report:
                    print(report, file=reportfh)
            outfile.write(segment)
