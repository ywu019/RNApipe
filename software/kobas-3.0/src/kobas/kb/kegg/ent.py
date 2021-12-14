#!/usr/bin/env python

import re

from kobas import utils

def parse(handle, abbr):
    genes, pathways, gene_pathways, ko_genes, \
    gene_entrez_gene_ids, gene_gis, gene_uniprotkb_acs, gene_ensembl_gene_ids, \
    ko_entrez_gene_ids, ko_gis, ko_uniprotkb_acs, ko_ensembl_gene_ids \
    = [], {}, [], {}, \
    {}, {}, {}, {}, \
    [], [], [], []

    is_pathway, is_ko, is_link, count = 0, 0, 0, 0
    for line in handle:
        if line.startswith('ENTRY'):
            gid = abbr + ':' + line.split()[1]
            genes.append([gid, ''])
        if line.startswith('NAME'):
            genes[-1][1] = line[:-1].split(None, 1)[1]

        if line.startswith('PATHWAY'):
            is_pathway = 1
        elif not line.startswith(' '):
            is_pathway = 0

        if line.startswith('ORTHOLOGY'):
            is_ko = 1
        elif not line.startswith(' '):
            is_ko = 0

        if line.startswith('DBLINKS'):
            is_link = 1
        elif not line.startswith(' '):
            is_link = 0

        if is_pathway:
            line = re.sub('PATHWAY', '', line)
            paid, pname = line[:-1].split(None, 1)
            if not pathways.has_key(pname):
                count += 1
                pathways[pname] = [count * 10 + 1, 'K', paid]
            gene_pathways.append((gid, pathways[pname][0]))

        if is_ko:
            line = re.sub('ORTHOLOGY', '', line)
            koid = line.split()[0]
            ko_genes.setdefault(koid, []).append(gid)

        if is_link:
            linkid = line.split()
            if 'NCBI-GeneID: ' in line:
                pos=linkid.index('NCBI-GeneID:')
                gene_entrez_gene_ids.setdefault(gid, []).extend(linkid[pos+1:])
            if 'NCBI-GI: ' in line:
                pos=linkid.index('NCBI-GI:')
                gene_gis.setdefault(gid, []).extend(linkid[pos+1:])
            if 'UniProt: ' in line:
                pos=linkid.index('UniProt:')
                gene_uniprotkb_acs.setdefault(gid, []).extend(linkid[pos+1:])
            if 'Ensembl: ' in line:
                pos=linkid.index('Ensembl:')
                gene_ensembl_gene_ids.setdefault(gid, []).extend(linkid[pos+1:])

    ko_entrez_gene_ids = utils.combine_dicts(ko_genes, gene_entrez_gene_ids)
    ko_gis = utils.combine_dicts(ko_genes, gene_gis)
    ko_uniprotkb_acs = utils.combine_dicts(ko_genes, gene_uniprotkb_acs)
    ko_ensembl_gene_ids = utils.combine_dicts(ko_genes, gene_ensembl_gene_ids)

    genes = utils.list_to_tuple(genes)
    pathways = utils.dict_to_list(pathways)
    gene_pathways = utils.list_to_tuple(gene_pathways)
    ko_genes = utils.split_values(ko_genes)
    gene_entrez_gene_ids = utils.split_values(gene_entrez_gene_ids)
    gene_gis = utils.split_values(gene_gis)
    gene_uniprotkb_acs = utils.split_values(gene_uniprotkb_acs)
    gene_ensembl_gene_ids = utils.split_values(gene_ensembl_gene_ids)

    return (genes, pathways, gene_pathways, ko_genes, \
    gene_entrez_gene_ids, gene_gis, gene_uniprotkb_acs, gene_ensembl_gene_ids, \
    ko_entrez_gene_ids, ko_gis, ko_uniprotkb_acs, ko_ensembl_gene_ids)

if __name__ == '__main__':
    import sys
    from pprint import pprint

    genes, pathways, gene_pathways, ko_genes, \
    gene_entrez_gene_ids, gene_gis, gene_uniprotkb_acs, gene_ensembl_gene_ids, \
    ko_entrez_gene_ids, ko_gis, ko_uniprotkb_acs, ko_ensembl_gene_ids \
    = parse(open(sys.argv[1]), sys.argv[2])

    print len(genes), len(pathways), len(gene_pathways), len(ko_genes)
    print len(gene_entrez_gene_ids), len(gene_gis), len(gene_uniprotkb_acs), len(gene_ensembl_gene_ids)
    print len(ko_entrez_gene_ids), len(ko_gis), len(ko_uniprotkb_acs), len(ko_ensembl_gene_ids)

    g, p = set(), set()
    for gp in gene_pathways:
        g.add(gp[0])
        p.add(gp[1])

    print len(g), len(p)

    pprint(genes[:5])
    pprint(pathways[:5])
    pprint(gene_pathways[:5])
    pprint(ko_genes[:5])
    pprint(gene_entrez_gene_ids[:5])
    pprint(gene_gis[:5])
    pprint(gene_uniprotkb_acs[:5])
    pprint(gene_ensembl_gene_ids[:5])
    pprint(ko_entrez_gene_ids[:5])
    pprint(ko_gis[:5])
    pprint(ko_uniprotkb_acs[:5])
    pprint(ko_ensembl_gene_ids[:5])
