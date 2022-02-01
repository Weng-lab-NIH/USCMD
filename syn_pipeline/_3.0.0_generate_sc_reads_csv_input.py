#!/usr/bin/env python3

from random import randint
from Bio.Seq import Seq
from Bio import SeqIO
import argparse
import csv
from math import floor, ceil
from pprint import pprint
from single_cell_classes import Donor, Cell, SequencedArea, UMI, ReadSet, Read
import gzip
from pathlib import Path
import os
import pandas as pd

READ_LENGTH=101
MAX_SKEW = int(max(10, READ_LENGTH/2 - 20))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('donor_exome', type=str)
    parser.add_argument('donor_info', type=str)
    parser.add_argument('donor_single_cell', type=str)
    parser.add_argument('barcode_whitelist', type=str)
    parser.add_argument('umi_whitelist', type=str)
    parser.add_argument('outdir', type=str)
    args = parser.parse_args()
    print(args)

    donor_info = pd.read_csv(args.donor_info)
    donor_single_cell = pd.read_csv(args.donor_single_cell)
    # print(args.donor_single_cell)
    # print(donor_single_cell)

    bc_whitelist_obj = open(args.barcode_whitelist)
    umi_whitelist_obj = open(args.umi_whitelist)
    donor_exome = list(SeqIO.parse(open(args.donor_exome), 'fasta'))
    donor_obj = generate_donor_object(donor_info,
        donor_single_cell, 
        bc_whitelist_obj, 
        umi_whitelist_obj,
        donor_exome)
    write_reads(donor_obj, args.outdir)

def generate_donor_object(donor_csv,
    single_cell_csv,
        bc_whitelist_obj, 
        umi_whitelist_obj,
        donor_exome):
    # print('generate_donor_object')
    donor_id = donor_csv.at[0,'donor_id']
    age = donor_csv.at[0, 'age']
    sex = donor_csv.at[0, 'sex']
    if (single_cell_csv['donor_id']!=donor_id).any():
        print("the donor_id column in donor_single_cell should match the" + 
              " donor_id" + 
              " in the first row of donor_info")
        exit(1)
    cell_objs = []
    # print(single_cell_csv)
    for cell_id in single_cell_csv['cell_id'].unique():
        print("cell_id:", cell_id)
        cell_df = single_cell_csv[single_cell_csv['cell_id'] == cell_id]
        cell_barcode = get_cell_barcode(bc_whitelist_obj)
        sequenced_area_objs = []
        for area_of_interest_id in cell_df['area_of_interest_id'].unique():
            AoI_df = cell_df[cell_df['area_of_interest_id'] == area_of_interest_id]
            AoI_df = AoI_df.reset_index(drop=True)
            sequenced_area_obj = generate_seq_area_obj(AoI_df, 
                cell_barcode, umi_whitelist_obj, donor_exome)
            sequenced_area_objs.append(sequenced_area_obj)
        cell_obj = Cell(
            cell_barcode = cell_barcode,
            sequenced_areas = sequenced_area_objs
            )
        cell_objs.append(cell_obj)
    donor_obj = Donor(
        donor_id = donor_id,
        age = age,
        sex = sex,
        cells = cell_objs
        )

    return donor_obj
    
def generate_seq_area_obj(csv_seq_area, cell_barcode, 
    umi_whitelist_obj, donor_exome):
    # print(csv_seq_area)
    umi_objs = []
    for field in ['exome_ref_line', 'position_of_interest']:
        if not len(csv_seq_area[field].unique()) == 1:
            print(csv_seq_area)
            print("this area of interest should have the same field in" + 
                "each row.")
            eixt(1)
    exome_ref_line = csv_seq_area.at[0, 'exome_ref_line']
    position_of_interest = csv_seq_area.at[0, 'position_of_interest']
    exome_subset = donor_exome[exome_ref_line].seq
    start_point = floor(position_of_interest - MAX_SKEW - READ_LENGTH/2)
    end_point = ceil(position_of_interest + MAX_SKEW + READ_LENGTH/2)
    exome_subset = list(exome_subset)[start_point : end_point]

    mut_loc_in_neighborhood = ceil(MAX_SKEW + READ_LENGTH/2)
    original_nt = exome_subset[mut_loc_in_neighborhood]
    for umi_num in csv_seq_area['UMI_id'].unique():
        umi_id = get_umi_id(umi_whitelist_obj)

        read_set_objs = []
        umi_df = csv_seq_area[csv_seq_area['UMI_id'] == umi_num]
        for read_set_id in umi_df['read_set_id'].unique():
            read_set_df = umi_df[umi_df['read_set_id'] == read_set_id].reset_index(drop=True)
            if read_set_df.shape[0] !=1:
                print(read_set_df)
                print("above df should have exactly 1 row")
                exit(1)

            print(f"mutation from {original_nt} to " + 
                f"{read_set_df.at[0, 'set_nucleotide_to']}")
            # print(read_set_df)
            exome_subset[mut_loc_in_neighborhood] = read_set_df.at[0,"set_nucleotide_to"]
            read_objs = []
            for read_idx in range(read_set_df.at[0,'num_reads']):
                read_start_point = randint(0, MAX_SKEW*2)
                read_end_point = read_start_point + READ_LENGTH
                # print(type(exome_subset))
                read_str = "".join(exome_subset[read_start_point : read_end_point])
                read_obj = Read(read_str)
                read_objs.append(read_obj)
            read_set_obj = ReadSet(reads = read_objs)
            read_set_objs.append(read_set_obj)
        umi_obj = UMI(
            UMI_id = umi_id,
            read_sets = read_set_objs
            )
        umi_objs.append(umi_obj)
    return SequencedArea(umi_objs)

def get_cell_barcode(bc_whitelist_obj):
    bc = bc_whitelist_obj.readline().strip()
    if len(bc) == 16:
        return bc
    print("ERROR: bc from bc_whitelist_obj is not 16 characters long.")
    print(bc, len(bc))
    exit(1)

def get_umi_id(umi_id_whitelist_obj):
    umi_id = umi_id_whitelist_obj.readline().strip()
    if len(umi_id) == 10:
        return umi_id
    print("ERROR: umi_id from umi_id_whitelist_obj is not 10 characters long.")
    print(umi_id, len(umi_id))
    exit(1)

def write_reads(donor_obj, outdir):
    used_barcode_path = f"{outdir}/used_barcodes.txt"
    if os.path.exists(used_barcode_path):
        os.remove(used_barcode_path)

    for cell in donor_obj.cells:
        barcode_w_dash = f"{cell.cell_barcode}-1"
        print("barcode", barcode_w_dash)
        with open(used_barcode_path, "a+") as f:
            f.write(f"{barcode_w_dash}\n")
        r1_records = []
        r2_records = []
        for sequenced_area in cell.sequenced_areas:
            for umi in sequenced_area.umis:
                print("umi id", umi.UMI_id)
                i = 0
                for read_set in umi.read_sets: 
                    for read in read_set.reads:
                        desc_str = ""
                        id_str = f"J00:sc_{barcode_w_dash}_{umi.UMI_id}"
                        quality_score1 = {"solexa_quality": [40 for character in read.sequence]}
                        r1_seq = Seq(read.sequence)
                        r1_record = SeqIO.SeqRecord(
                            r1_seq, 
                            id=id_str,
                            letter_annotations=quality_score1,
                            description = desc_str
                            )
                        r1_records.append(r1_record)

                        quality_score2 = {"solexa_quality": [40 for character in range(26)]}
                        r2_seq = Seq(cell.cell_barcode + umi.UMI_id)
                        r2_record = SeqIO.SeqRecord(
                            r2_seq, 
                            id=id_str,
                            letter_annotations=quality_score2,
                            description = desc_str
                            )
                        r2_records.append(r2_record)
                        i+=1

        r1_out = str(Path(outdir) / f"{barcode_w_dash}_r1.fastq.gz")
        r1_out = gzip.open(r1_out, "wt")
        r2_out = str(Path(outdir) / f"{barcode_w_dash}_r2.fastq.gz")
        r2_out = gzip.open(r2_out, "wt")
        SeqIO.write(r1_records, handle=r1_out, format='fastq')
        SeqIO.write(r2_records, handle=r2_out, format='fastq')

if __name__ == "__main__":
    main()
