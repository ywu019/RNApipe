#!/usr/bin/env python

from kobas import utils

def parse(map_handle, term_handle, kobasdb):
    diseases, gene_diseases = {}, set()

    count = 0
    for line in term_handle:
        if line.startswith('    <li class="'):
            line = line.replace("&quot;", "")
            line = line.replace("&#39;", "'")

            count += 1
            diid, dname = line.split('objects/')[1].split('</a>')[0].split('/">')
            diseases[dname] = [count * 10 + 3, 'f', diid]

    for line in map_handle:
        dname, entrez_gene_id = line.split('\t')[0:2]
        dname = dname.replace('"', '')
        gids = kobasdb.gids_from_entrez_gene_id(entrez_gene_id)
        for gid in gids:
            gene_diseases.add((gid[0], diseases[dname][0]))

    dids = set()
    for gene_disease in gene_diseases:
        dids.add(gene_disease[1])

    for disease in diseases.items():
        if disease[1][0] not in dids:
            diseases.pop(disease[0])

    diseases = utils.dict_to_list(diseases)
    gene_diseases = list(gene_diseases)

    return diseases, gene_diseases

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config, dbutils

    diseases, gene_diseases = parse(open(sys.argv[1]), open(sys.argv[2]), \
        dbutils.KOBASDB(config.getrc()['kobasdb'] + 'hsa.db'))

    print len(diseases), len(gene_diseases)

    g, d = set(), set()
    for gd in gene_diseases:
        g.add(gd[0])
        d.add(gd[1])

    print len(g), len(d)

    pprint(diseases[:5])
    pprint(gene_diseases[:5])
