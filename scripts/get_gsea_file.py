import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to get GSEA prepared files.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", metavar="--input", required=True, type=str, help="input file about geneid_tsv")
    parser.add_argument("-s", metavar="--Samples", required=True, type=str, help="samples_tsv")
    parser.add_argument("--treat", required=True, type=str, help="treat:['KO']")
    parser.add_argument("--control", required=True, type=str, help="control:['WT']")
    parser.add_argument("-e", metavar="--entrezid", required=True, type=str, help="EntrezID_list")
    parser.add_argument("-g", metavar="--gct", required=True, type=str, help="output_gct_file(.txt)")
    parser.add_argument("-c", metavar="--cls", required=True, type=str, help="output_cls_file")
    args = parser.parse_args()
    return args

def get_gct(geneid_tsv, entrezid, gct_tsv):
    entrez_id = pd.read_csv(entrezid, sep="\t", header=None)
    entrez_id.columns = ["Geneid","EntrezID"]
    gene_id_tsv = pd.read_csv(geneid_tsv, sep="\t")
    merge_data = pd.merge(gene_id_tsv,entrez_id,how="inner",on="Geneid")
    df = merge_data.reindex(columns=["EntrezID","Geneid"]+list(gene_id_tsv.columns[1:]))
    df=df.rename(columns={"Geneid":"description","EntrezID":"Geneid"})
    df.to_csv(gct_tsv, sep="\t", index=False)


def get_cls(Samples, cls_file, treat, control):
    samples = pd.read_csv(Samples, sep="\t")
    line1 = str(len(samples)) + " 2 1"
    line2 = "#" + " " + str(treat) + " " + str(control)
    line3 = " ".join(samples["condition"])
    with open(cls_file, 'w') as file_out:
        file_out.writelines(line1 + "\n")
        file_out.writelines(line2 + "\n")
        file_out.writelines(line3)

def main():
    args = parse_args()
    get_gct(args.i, args.e, args.g)
    get_cls(args.s, args.c, args.treat, args.control)

if __name__ == "__main__":
    main()

