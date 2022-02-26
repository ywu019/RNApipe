import numpy as np
import pandas as pd
import argparse
import os
import re

def parse_args():
    parser = argparse.ArgumentParser(description = "A feature tpm.")
    parser.add_argument('-i', '--countsfile', action='store', dest='countsfile')
    parser.add_argument('-o1', '--outfile_count', action='store',dest='outfile_count')
    parser.add_argument('-o2', '--outfile_tpm', action='store',dest='outfile_tpm')
    args = parser.parse_args()
    return args

def counts_to_tpm(count,sizes):
    rate = np.log(count).subtract(np.log(sizes))
    denom = np.log(np.sum(np.exp(rate)))
    tpm = np.exp(rate - denom + np.log(1e6))
    return tpm

def featurecounts2tpm(countsfile, outfile_count, outfile_tpm):
    counts = pd.read_csv(countsfile,sep="\t")
    counts = counts.set_index("Geneid")
    counts = counts.drop(columns=['Chr','Start','End','Strand'])
    lengths=counts['Length']
    counts = counts.drop(columns=['Length'])
    counts.columns = [re.sub(".sort.*","",os.path.basename(col)) for col in counts.columns]
    counts.to_csv(outfile_count,sep="\t",index=True,header=True)
    tpm = counts.apply(lambda x: counts_to_tpm(x,lengths),axis=0)
    tpm.to_csv(outfile_tpm,sep="\t",index=True,header=True)

def main():
    args = parse_args()
    featurecounts2tpm(args.countsfile, args.outfile_count, args.outfile_tpm)

if __name__ == '__main__':
    main()

