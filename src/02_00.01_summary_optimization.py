#!/bin/python3

import sys
import os
import re
import pysam

def parse_mia_file(file):
    """
    Parse the output of MIA (Mapping-Iterative-Assembler) and return the ID,
    staring and end positions of each mapped read.
    NOTE: MIA only stores mapped reads.
    """
    MIA_ID = []
    MIA_START = []
    
    with open(file, 'r', errors='replace') as f:
        lines = f.readlines()[10:]
        for line in lines:
            if line[0:2] == "ID":
                MIA_ID.append(line.strip()[3:].rstrip("_b").rstrip("_t"))
            elif line[0:5] == "START":
                MIA_START.append(int(line.strip()[6:]))
    
    return MIA_ID, MIA_START


def process_bam_file(file, MIA_ID, MIA_START):
    TP = 0
    TN = 0
    FP = 0
    FN = 0
    
    samfile = pysam.AlignmentFile(file, "rb", index_filename=file+".bai")
    for read in samfile.fetch(until_eof=True):
        # Get info of each read
        ID_i = read.query_name
        FLAG_i = read.flag
        START_i = read.pos
        MAPQ = read.mapq
        # if ID_i is in MIA_ID means that the read i was mapped with MIA
        if ID_i in MIA_ID:
            # if FLAG_i is unmapped or MAPQ is low, then it is a false negative
            if FLAG_i == 4 or MAPQ < 20:
                FN += 1
            # else, could be a match or a mistake
            else:
                idx = MIA_ID.index(ID_i)
                # consider true positive only if start is the same position or
                # within 35 nt
                if abs(START_i - MIA_START[idx]) <= 35:
                    TP += 1
                # else, it is a false positive
                else:
                    FP += 1
        # else, ID_i was not mapped with MIA
        else:
            # if read ID_i is unmapped or low MAPQ, it is a true negative
            if FLAG_i == 4 or MAPQ < 20:
                TN += 1
            # else, the read was mapped somewere; hence it is a false positive
            else:
                FP += 1
    
    return TP, FP, TN, FN

def get_metrics(TP, FP, TN, FN):
    acc = (TP + TN) / (TP + TN + FP + FN)
    sens = TP / (TP + FN)
    spec = TN / (TN + FP)
    ba_acc = (sens + spec) / 2
    f1 = (2*TP)/(2*TP + FP + FN)
    return [TP, FP, TN, FN, acc, sens, spec, ba_acc, f1]

if __name__ == '__main__':
    output = open("results/00_hDNA_mapper_opt/metrics.tsv", "w+")
    output.write("file\tTP\tFP\tTN\tFN\tacc\tsens\tspec\tba_acc\tf1\n")
    # parse the MIA file
    MIA_ID, MIA_START = parse_mia_file("results/00_hDNA_mapper_opt/00_MIA/MIA_assembly.4")
    # get the bam file list
    file_list = os.listdir("results/00_hDNA_mapper_opt")
    ptrn = re.compile('^.*.bam$')
    bam_list = [ s for s in file_list if ptrn.match(s)]
    # parse bam files and get metrics
    for file in bam_list:
        TP, FP, TN, FN = process_bam_file("results/00_hDNA_mapper_opt/"+file,
                                          MIA_ID, MIA_START)
        metrics = get_metrics(TP, FP, TN, FN)
        output.write(file+"\t"+"\t".join([str(m) for m in metrics])+"\n")
    output.close()
