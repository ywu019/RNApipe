#!/usr/bin/env python

import re

from kobas import utils

def parse(handle):
    diseases, gene_diseases = {}, []

    count = 0
    for line in handle:
        if line.startswith('A'):
            adname = re.sub('.*<b>(.*)</b>.*', lambda mo: mo.group(1), line).strip()
            count += 1
            diseases[adname] = [count * 10 + 2, 'k', '']
        if line.startswith('B'):
            bdname = line[:-1].split(None, 1)[1]
            if not diseases.has_key(bdname):
                count += 1
                diseases[bdname] = [count * 10 + 2, 'k', '']
        if line.startswith('C'):
            cdiid, cdname = line.split('[')[0].strip().split(None, 2)[1:]
            if '-' in cdiid:
                cdiid = ''
            if not diseases.has_key(cdname):
                count += 1
                diseases[cdname] = [count * 10 + 2, 'k', cdiid]
        if line.startswith('D'):
            if '[HSA:' in line:
                gidroots = re.sub('.*\[HSA:(.*?)\].*', lambda mo: mo.group(1), line).split()
                for gidroot in gidroots:
                    gid = 'hsa:' + gidroot
                    gene_diseases.append((gid, diseases[adname][0]))
                    gene_diseases.append((gid, diseases[bdname][0]))
                    gene_diseases.append((gid, diseases[cdname][0]))

    dids = set()
    for gene_disease in gene_diseases:
        dids.add(gene_disease[1])

    for disease in diseases.items():
        if disease[1][0] not in dids:
            diseases.pop(disease[0])

    diseases = utils.dict_to_list(diseases)
    gene_diseases = utils.list_to_tuple(gene_diseases)

    return (diseases, gene_diseases)

if __name__ == '__main__':
    import sys
    from pprint import pprint

    diseases, gene_diseases = parse(open(sys.argv[1]))

    print len(diseases), len(gene_diseases)

    g, d = set(), set()
    for gd in gene_diseases:
        g.add(gd[0])
        d.add(gd[1])

    print len(g), len(d)

    pprint(diseases[:5])
    pprint(gene_diseases[:5])
