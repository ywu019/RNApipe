#!/usr/bin/env python

from kobas import utils
from kobas.kb.nhgri import gene_info, correction

def parse(gwas_handle, gene_handle, kobasdb):
    diseases, gene_diseases = {}, set()

    entrez_gene_id_from_symbol = gene_info.parse(gene_handle)

    count = 0
    gwas_handle.readline()    ##skip the first line
    for line in gwas_handle:
        line = line[:-1].split('\t')
        if len(line) >= 14:
            dname, symbols = line[7].strip(), line[13].strip()
            dname = dname.replace("&beta;", "beta")
            if correction.correction.has_key(symbols):
                symbols = correction.correction[symbols]

            if dname:
                for symbol in symbols.split(','):
                    symbol = symbol.strip()
                    if entrez_gene_id_from_symbol.has_key(symbol):
                        entrez_gene_id = entrez_gene_id_from_symbol[symbol]
                        gids = kobasdb.gids_from_entrez_gene_id(entrez_gene_id)
                        for gid in gids:
                            if not diseases.has_key(dname):
                                count += 1
                                diseases[dname] = [count * 10 + 5, 'N', '']

                            gene_diseases.add((gid[0], diseases[dname][0]))

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
