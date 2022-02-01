#!/usr/bin/env python3

import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sam_in', type=str)
    parser.add_argument('barcode', type=str)
    parser.add_argument('sam_out', type=str)
    args = parser.parse_args()

    print(args)

    if os.path.exists(args.sam_out):
        os.remove(args.sam_out)

    sam_in = open(args.sam_in, "r")
    sam_out = open(args.sam_out, "a+")

    for line_in in sam_in:
        line_list = list(line_in.strip())
        line_words = line_in.split()
        if line_list[0] == 'J':
            umi_barcode = line_words[0].split("_")[2]
            ub_tag = f"\tUB:Z:{umi_barcode}"
            cb_tag = f"\tCB:Z:{args.barcode}"
            line_list.extend(list(cb_tag))
            line_list.extend(list(ub_tag))
        line_list.append('\n')
        line_out = "".join(line_list)
        #print(f"from {line_in} to {line_out}")

        sam_out.write(line_out)

if __name__ == "__main__":
    main()
