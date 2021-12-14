#!/usr/bin/env python

from kobas import utils

def parse(handle, kobasdb):
    diseases, gene_diseases = {}, set()

    count = 0
    for line in handle:
        line = line[:-1].split('\t')
        if len(line) >= 23:
            association, broad_phenotype, disease_class, entrez_gene_id = line[1:4] + line[22:23]
            if association == 'Y':
                if broad_phenotype:
                    if '|' in broad_phenotype:
                        phenotypes = broad_phenotype.split('|')
                    elif ';' in broad_phenotype:
                        phenotypes = broad_phenotype.split(';')
                    else:
                        phenotypes = [broad_phenotype]
                else:
                    phenotypes = []

                if disease_class not in ['OTHER', 'UNKNOWN']:
                    phenotypes.append(disease_class)

                gids = kobasdb.gids_from_entrez_gene_id(entrez_gene_id)

                for gid in gids:
                    for phenotype in phenotypes:
                        dname = phenotype.strip().upper()
                        if dname.endswith('.'):
                            dname = dname[:-1]
                        if dname:
                            if not diseases.has_key(dname):
                                count += 1
                                diseases[dname] = [count * 10 + 4, 'g', '']

                            gene_diseases.add((gid[0], diseases[dname][0]))

    diseases = utils.dict_to_list(diseases)
    gene_diseases = list(gene_diseases)

    return diseases, gene_diseases

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config, dbutils

    diseases, gene_diseases = parse(open(sys.argv[1]), \
        dbutils.KOBASDB(config.getrc()['kobasdb'] + 'hsa.db'))

    print len(diseases), len(gene_diseases)

    g, d = set(), set()
    for gd in gene_diseases:
        g.add(gd[0])
        d.add(gd[1])

    print len(g), len(d)

    pprint(diseases[:5])
    pprint(gene_diseases[:5])
