import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="This script is used to get gene length list.",
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", required=True, type=str, help="in_file about gtf file")
    parser.add_argument("-o", required=True, type=str, help="output file about gene length list")
    args = parser.parse_args()
    return args

def get_longest_trans(in_file, out_file):
    gtf = pd.read_csv(in_file, sep='\t',
                      names=["Chr", "database", "type", "start", "end", "score", "strand", "phase", "description"])
    gtf = gtf[~gtf["Chr"].str.contains("#")].reset_index(drop=True)
    gtf_exon = gtf[gtf["type"].str.contains('exon')]
    start = list(gtf_exon["start"])
    end = list(gtf_exon["end"])
    des = list(gtf_exon["description"])
    length_exon = []
    for i in range(len(start)):
        length = int(end[i]) - int(start[i]) + 1
        length_exon.append(length)
    length_list = []
    longest_trans = {}
    with open(out_file, "w") as file_out:
        for index, value in enumerate(des):
            name_Trans = value.split(';')[1].split(' ')[-1].strip('"')
            name_Genes = value.split(';')[0].split(' ')[-1].strip('"')
            new = length_exon[index]
            if not longest_trans.get(name_Genes):
                longest_trans[name_Genes] = {}
                longest_trans[name_Genes][name_Trans] = int(new)
            elif longest_trans.get(name_Genes) and longest_trans[name_Genes].get(name_Trans):
                longest_trans[name_Genes][name_Trans] = longest_trans[name_Genes][name_Trans] + int(new)
            elif longest_trans.get(name_Genes) and not longest_trans[name_Genes].get(name_Trans):
                longest_trans[name_Genes][name_Trans] = int(new)
        for gene in longest_trans.keys():
            for trans in longest_trans[gene].keys():
                length_list.append((trans, longest_trans[gene][trans]))
            length_list.sort(key=lambda x: x[1])
            longest = str(length_list[-1][1])
            longest = gene + "\t" + longest
            file_out.writelines(longest + "\n")
            length_list = []



def main():
    args = parse_args()
    get_longest_trans(args.i, args.o)

if __name__ == "__main__":
    main()

