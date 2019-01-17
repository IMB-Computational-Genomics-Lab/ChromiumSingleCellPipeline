#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
from collections import OrderedDict
import os
import sys
import argparse

def writeGTF(output_lines, output_filename):
    if len(output_lines) > 1:
        output_lines = ("\n").join(output_lines)
    else:
        output_lines = output_lines[0]

    with open(output_filename, "w") as output:
        output.write(output_lines)
    print ("GTF reference generated!")

def formatGTF(fasta_dict):
    # This contains three sets of entries for a GTF file
    gtf_template = """{genome}\thavana\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{id}"; gene_name "{id}"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n{genome}\thavana\ttranscript\t{start}\t{end}\t.\t+\t.\tgene_id "{id}"; transcript_id "{id}"; gene_name "{id}"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "{id}"; transcript_source "havana";\n{genome}\thavana\texon\t{start}\t{end}\t.\t+\t.\tgene_id "{id}"; transcript_id "{id}"; exon_number "1"; gene_name "{id}"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "{id}"; transcript_source "havana"; exon_id "{id}_exon";"""

    output_lines = []

    for gene, seq in fasta_dict.items():
        variable_dict = {"genome": gene, "start": 1, "end": len(seq), "id": gene}
        formatted_string = gtf_template.format(**variable_dict)
        output_lines.append(formatted_string)

    return output_lines

def parseFASTA(input_filename):
    # Parse FASTA to collect sequences and genes
    with open(input_filename, "r") as input:
        # Collect the reads in a file
        read_dict = OrderedDict()
        chrom_name = ""

        # Iterate through each line in the fasta file
        for line in input:
            line = line.strip()
            # If this is a fasta header, grab the information
            if line.startswith(">"):
                chrom_name = line.replace(">", "")

                if chrom_name not in read_dict.keys():
                    read_dict[chrom_name] = []
            else:
                current_reads = read_dict[chrom_name]
                current_reads.append(line)

    # Iterate through collected genes to concatenate FASTA sequences
    output_dict = OrderedDict()
    for gene, sequence in read_dict.items():
        sequence = ("").join(sequence)
        output_dict[gene] = sequence.upper()
    return output_dict

def parseArgs():
    parser = argparse.ArgumentParser(prog = "nucFasta2GTF.py", description = "Generates GTF files from a FASTA sequence.")
    parser.add_argument("-i", "--input", type = str, help = "A nucleotide sequence in FASTA format.", required = True)
    parser.add_argument("-o", "--output", type = str, help = "Output GTF filename.", required = True)
    args = parser.parse_args()
    return args.input, args.output

if __name__ == "__main__":
    input_filename, output_filename = parseArgs()

    # Parse the FASTA file to collect transcript name and dimensions
    transcript_data = parseFASTA(input_filename)

    # Use FASTA data to generate GTF
    output_data = formatGTF(transcript_data)
    writeGTF(output_data, output_filename)
