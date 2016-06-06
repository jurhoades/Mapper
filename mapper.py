#!/usr/bin/env python

# Imports
from __future__ import print_function
import argparse
import itertools
import math
import os
import sys

# Functions
def parse_arguments():
    """Parse command line arguments.

    Parse all command line arguments and return parameters.
    """
    parser = argparse.ArgumentParser(description="Map genomic coordinate to specific RefSeq ID CDS and AA position.")
    parser.add_argument("RefSeq_ID", type=str, help="RefSeq ID (ex NM_146145)")
    parser.add_argument("coordinate", type=int, help="genomic coordinate")
    parser.add_argument("filename", type=check_file, 
                        help="path to refFlat file used for mapping")
    parser.add_argument("-v", "--verbose", help="increase text output", 
                        action="store_true")

    return parser.parse_args()


def check_file(filename):
    if not os.path.exists(filename):
        msg = "Invalid file path"
        raise argparse.ArgumentTypeError(msg)

    return os.path.abspath(filename)


# Classes
class Mapper(object):
    
    def __init__(self, filename):
        self.__refflat_idx = self.__index_refflat(filename)

    
    def map(self, coord, refseq_id):
        try:
            #get gene annotation
            gene_anno = self.__get_annotation(refseq_id)

            #create gene structure
            gene_structure = self.__create_structure(gene_anno)
    
            #map genomic coordinate to CDS and AA position
            cds_pos, aa_pos = self.__find_position(gene_structure, coord)

            return cds_pos, aa_pos

        except Exception:
            return None, None


    def __index_refflat(self, refflat_file):
        with open(refflat_file) as read_file:
            refflat_idx = {}
            for line in read_file:
                line = line.strip().split()
                refseq_id = line[1]
                if refseq_id[:2] == "NM":
                    refflat_idx[refseq_id] = line

        return refflat_idx


    def __get_annotation(self, refseq_id):
        return self.__refflat_idx.get(refseq_id)


    def __create_structure(self, anno):
        """Create mapping structure from gene annotation.

        Use refFlat gene annotation to create a list of genomic coordinates
        for a particular coding sequence.
        """
        exon_starts = anno[9][:-1].split(",")
        exon_ends = anno[10][:-1].split(",")

        #trim exons not in CDS
        while anno[6] > exon_ends[0]:
            del exon_starts[0]
            del exon_ends[0]
        while anno[7] < exon_starts[-1]:
            del exon_starts[-1]
            del exon_ends[-1]

        exon_starts[0] = anno[6]
        exon_ends[-1] = anno[7]

        coding_seq_ranges = [[exon_starts[i], exon_ends[i]] for i in range(len(exon_starts))]
        coding_seq = [list(range(int(rng[0])+1,int(rng[1])+1)) for rng in coding_seq_ranges]
        coding_seq = list(itertools.chain(*coding_seq))  ##combine exon list of lists

        if len(coding_seq) % 3 != 0:
            raise Exception

        structure = {"coding_seq": coding_seq, "strand": anno[3]}

        return structure


    def __find_position(self, structure, genome_loc):
        """Find CDS and AA position from gene structure.

        Uses gene structure to determine CDS and amino acid position of
        genomic location.
        """
        strand = structure["strand"]
        coding_seq = structure["coding_seq"]

        idx = coding_seq.index(genome_loc)

        if strand == "+":
            cds_pos = idx + 1
            aa_pos = (idx // 3) + 1
        elif strand == "-":
            cds_pos = len(coding_seq) - idx
            aa_pos = int(math.ceil(cds_pos / 3))

        return cds_pos, aa_pos        

        
# Main
def main():
    args = parse_arguments()

    mapper = Mapper(args.filename)

    cds_pos, aa_pos = mapper.map(args.coordinate, args.RefSeq_ID)

    if args.verbose:
        print("RefSeq ID: {0}\nGenomic Coordinate: {1}\nCDS Position: {2}\nAmino Acid Position: {3}".format(
              args.RefSeq_ID, args.coordinate, cds_pos, aa_pos))
    else:
        print(cds_pos, aa_pos)


if __name__ == "__main__":
    main()
