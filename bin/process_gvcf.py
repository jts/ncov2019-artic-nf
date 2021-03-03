#!/usr/bin/env python

import itertools
import argparse
import pysam
import sys
import os
from collections import defaultdict

# from https://www.geeksforgeeks.org/python-make-a-list-of-intervals-with-sequential-numbers/
# via artic-mask
def intervals_extract(iterable):
    iterable = sorted(set(iterable))
    for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
        group = list(group)
        yield [group[0][1], group[-1][1]]

def write_depth_mask(out_filename, contig_depths, min_coverage):
    maskfh = open(out_filename, 'w')
    for contig_name, depths in contig_depths.items():
        # from artic-mask, create list of positions that fail the depth check
        mask_vector = []
        for pos, depth in enumerate(depths):
            if depth < min_coverage:
                mask_vector.append(pos)

        # get the intervals from the mask_vector
        intervals = list(intervals_extract(mask_vector))

        for i in intervals:
            maskfh.write("%s\t%s\t%s\n" % (contig_name, i[0]+1, i[1]+1))
    maskfh.close()

def main():

    description = 'Process a .gvcf file to create a file of consensus variants, low-frequency variants and a coverage mask'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-m', '--mask-output', required=True,
            help=f"The output file name for the coverage mask\n")
    
    parser.add_argument('-v', '--variants-output', required=True,
            help=f"The output file name for variants (non-reference gVCF records)\n")

    parser.add_argument('-d', '--min-depth', type=int, default=10,
            help=f"Mask reference positions with depth less than this threshold")
    
    parser.add_argument('-l', '--lower-ambiguity-frequency', type=float, default=0.25,
            help=f"Variants with frequency less than -l will be discarded")
    
    parser.add_argument('-u', '--upper-ambiguity-frequency', type=float, default=0.75,
            help=f"Substitution variants with frequency less than -u will be encoded with IUPAC ambiguity codes")
    
    parser.add_argument('file', action='store', nargs=1)
    
    args = parser.parse_args()
    vcf = pysam.VariantFile(open(args.file[0],'r'))

    # Initalize depth mask to all zeros for all contigs
    contig_depth = defaultdict(list)
    for r in vcf.header.records:
        if r.type == "CONTIG":
            contig_depth[r['ID']] = [0] * int(r['length'])

    out_header = vcf.header
    out_header.info.add("VAF", number="A", type='Float', description="Variant allele fraction, called from observed reference/alt reads")
    out_header.info.add("ConsensusTag", number=1, type='String', description="The type of base to be included in the consensus sequence (IUPAC or Fixed)")
    variants_out = pysam.VariantFile(args.variants_output,'w',header=out_header)

    for record in vcf:

        is_gvcf_ref = record.alts[0] == "<*>"
        #print(is_gvcf_ref, record.alts[0])

        # set depth for this part of the genome
        # this works for both gVCF blocks and regular variants
        # because pos/stop are set appropriately
        v_start = record.pos
        v_end = record.stop
        depth = record.info["DP"]

        # disallow gvcf records that are longer than a single base
        assert(not is_gvcf_ref or v_start == v_end)

        for i in range(v_start, v_end + 1):
            assert(i > 0)
            # VCF coordinates are 1-based, we record the depth vector as 0-based
            # to be consistent with artic-mask
            contig_depth[record.chrom][i - 1] = depth

        # do nothing else with ref records
        if is_gvcf_ref:
            continue

        # calculate VAF for each allele
        vafs = list()
        total_depth = float(record.info["DP"])
        sum_vaf = 0.0
        max_vaf = 0.0
        max_idx = -1
        for i in range(0, len(record.alts)):
            alt_reads = int(record.info["AO"][i])
            vaf = float(alt_reads) / float(record.info["DP"])
            if vaf > max_vaf:
                max_vaf = vaf
                max_idx = i

            vafs.append(vaf)
            sum_vaf += vaf

        # check the assumption that freebayes writes the alleles in VAF order
        assert(max_idx == 0)

        # discard low frequency
        if vafs[max_idx] < args.lower_ambiguity_frequency or depth < args.min_depth:
            continue
    
        is_indel = len(record.ref) != len(record.alts[0])
        record.info["VAF"] = vafs
        consensus_tag = "None"

        # Write a tag describing what to do with the variant
        # For indels we use the summed VAF across all alleles to handle multi-allelic indels (e.g. CTTT ->  CTT, CT)
        if vafs[max_idx] > args.upper_ambiguity_frequency or (is_indel and sum_vaf >= 0.5):
            # always apply these to the consensus
            consensus_tag = "fixed"
        elif is_indel:
            # we can't represent ambiguous indels so they get this tag and are not applied to the consensus
            consensus_tag = "low_frequency_indel"
        else:
            # record ambiguous SNPs in the consensus sequence with IUPAC codes
            consensus_tag = "ambiguous"
        record.info["ConsensusTag"] = consensus_tag
        variants_out.write(record)

    write_depth_mask(args.mask_output, contig_depth, args.min_depth)

if __name__ == "__main__":
    main()
