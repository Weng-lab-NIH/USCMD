#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('donor_info', type=str)
    parser.add_argument('SNP_info', type=str)
    parser.add_argument('ref_path', type=str)
    parser.add_argument('outdir', type=str)
    args = parser.parse_args()

    donor_info = pd.read_csv(args.donor_info)
    SNP_info = pd.read_csv(args.SNP_info)
    ref_path = args.ref_path

    # line_number, "nucleotide number", "new nucleotide", "num_reads"
    # remeber that there are header lines in here!    
    donor_id = donor_info.at[0,'donor_id']
    if (SNP_info['donor_id']!=donor_id).any():
        print("the donor_id column in SNP_info should match the donor_id " + 
              " in the first row of donor_info")

    for allele_idx in [1,2]:
        allele_specific = SNP_info[SNP_info['allele_idx'] == allele_idx]
        ref = SeqIO.parse(open(ref_path),'fasta')
        new_exome = list(ref)
        out_path = str(
            Path(args.outdir) / f"{donor_id}_allele_{allele_idx}.fasta")

        for index,row in allele_specific.iterrows():
            temp_seq = new_exome[row['exome_ref_line']].seq
            position = row['position']
            new_base = row['set_nucleotide_to']
            temp_list = list(temp_seq)

            temp_subset = temp_list[position-100 : position+100]
            if 'N' in temp_subset:
                print("it contains N!")
            print(f"from {temp_list[position]} to {new_base}")
            if temp_list[position] == new_base:
                print("THIS ISN'T REALLY A SNP")
            temp_list[position] = new_base
            temp_seq = Seq("".join(temp_list))
            new_exome[row['exome_ref_line']].seq = temp_seq
            print("--------------")

        with open(out_path, 'w') as f:
            for line in new_exome:
                #print("type(line)", type(line))
                f.write(line.format("fasta"))

if __name__ == "__main__":
    main()
