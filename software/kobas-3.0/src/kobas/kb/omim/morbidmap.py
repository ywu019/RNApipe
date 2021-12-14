#!/usr/bin/env python

from kobas import utils
from kobas.kb.omim import correction

def parse(omim_handle, gene_handle, kobasdb):
    diseases, gene_diseases = {}, set()

##Disorder|GeneSymbols|GeneMIM|
    disorder_genes = {}
    for line in omim_handle:
        line = line.strip()
        if correction.correction.has_key(line):
            line = correction.correction[line]
        disorder, tmp, gene = line.split('|')[:3]

        disorder = disorder.rsplit(' ', 1)[0]
        if disorder[-6:].isdigit():    ##remove disorders without OMIM entry ID
            dname, diid = disorder.rsplit(',', 1)[0], disorder[-6:]
            if (dname[0] == '[' and dname[-1] == ']') or (dname[0] == '{' and dname[-1] == '}'):
                dname = dname[1:-1]
            if dname[0] == '?':
                dname = dname[1:]
            diseases.setdefault(diid, set()).add(dname)
            disorder_genes.setdefault(diid, set()).add(gene)

##tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\tdbXrefs\t
    gene_entrez_gene_ids = {}
    for line in gene_handle:
        entrez_gene_id, tmp1, tmp2, tmp3, dbxrefs = line.split('\t')[1:6]
        if 'MIM' in dbxrefs:
            mim_i = dbxrefs.find('MIM')
            gene = dbxrefs[mim_i + 4: mim_i + 10]
            gene_entrez_gene_ids[gene] = entrez_gene_id

    count = 0
    for item in diseases.items():
        count += 1
        dnames = list(item[1])
        diseases[item[0]] = [count * 10 + 1, 'o', dnames[0]]
        for dname in dnames[1:]:
            if len(dname) < len(diseases[item[0]][2]):
                diseases[item[0]][2] = dname

        if disorder_genes.has_key(item[0]):
            genes = disorder_genes[item[0]]
            for gene in genes:
                if gene_entrez_gene_ids.has_key(gene):
                    entrez_gene_id = gene_entrez_gene_ids[gene]
                    gids = kobasdb.gids_from_entrez_gene_id(entrez_gene_id)
                    for gid in gids:
                        gene_diseases.add((gid[0], diseases[item[0]][0]))

    dids = set()
    for gene_disease in gene_diseases:
        dids.add(gene_disease[1])

    for disease in diseases.items():
        if disease[1][0] not in dids:
            diseases.pop(disease[0])

    diseases = utils.dict_to_list_reorder(diseases)
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
