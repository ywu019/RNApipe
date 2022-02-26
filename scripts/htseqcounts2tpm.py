import numpy as np
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="htseq count to tpm.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', required=True, type=str, help="countsfile")
    parser.add_argument('-l', required=True, type=str, help="gene_length list")
    parser.add_argument('-o', required=True, type=str, help="outfile_tpm")
    args = parser.parse_args()
    return args


def counts_to_tpm(count, sizes):
    rate = np.log(count).subtract(np.log(sizes))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    return tpm


def htseq_counts2tpm(countfile, length_file, outfile_tpm):
    counts = pd.read_csv(countfile, sep="\t")
    lengths = pd.read_csv(length_file, sep="\t", header=None)
    lengths.columns = ["Geneid", "length"]
    count_tab = pd.merge(counts, lengths, on="Geneid", how="outer")
    count_tab = count_tab.set_index("Geneid")
    lengths = count_tab["length"]
    count_tab = count_tab.drop(columns=['length'])
    tpm = count_tab.apply(lambda x: counts_to_tpm(x, lengths), axis=0)
    tpm.to_csv(outfile_tpm, sep="\t", index=True, header=True)


def main():
    args = parse_args()
    htseq_counts2tpm(args.i, args.l, args.o)


if __name__ == "__main__":
    main()
