#!/usr/bin/env python

import re

def parse(gene_handle, pathway_handle, kobasdb):
    pathways, gene_pathways = set(), set()

    bcids = {}
    for line in gene_handle:
        if line.startswith('#'):
            pass
        elif line.startswith('UNIQUE-ID'):
            line = line[:-1].split('\t')
            bcid_index, uniprotkb_ac_index = line.index('UNIQUE-ID'), line.index('SWISS-PROT-ID')
        else:
            line = line[:-1].split('\t')
            bcid, uniprotkb_ac = line[bcid_index], line[uniprotkb_ac_index]
            if uniprotkb_ac:
                gids = kobasdb.gids_from_uniprotkb_ac(uniprotkb_ac)
                for gid in gids:
                    bcids.setdefault(bcid, set()).add(gid[0])

    count = 0
    for line in pathway_handle:
        if line.startswith('#'):
            pass
        elif line.startswith('UNIQUE-ID'):
            line = line[:-1].split('\t')
            paid_index, pname_index, gene_id_start = line.index('UNIQUE-ID'), line.index('NAME'), line.index('GENE-ID')
            gene_id_end = gene_id_start + line.count('GENE-ID')
        else:
            line = line[:-1].split('\t')
            paid, pname, gene_ids = line[paid_index], line[pname_index], line[gene_id_start:gene_id_end]
            count += 1
            pid = count * 10 + 5

            if '<' in pname:
                pname = re.sub('<i>(.+?)</i>', lambda mo: mo.group(1), pname)
                pname = re.sub('<I>(.+?)</I>', lambda mo: mo.group(1), pname)
                pname = re.sub('<sub>(.+?)</sub>', lambda mo: mo.group(1), pname)
                pname = re.sub('<SUB>(.+?)</SUB>', lambda mo: mo.group(1), pname)
                pname = re.sub('<sup>(.+?)</sup>', lambda mo: mo.group(1), pname)
            elif '&' in pname:
                pname = re.sub('&(.+?);', lambda mo: mo.group(1), pname)

            for gene_id in gene_ids:
                if bcids.has_key(gene_id):
                    pathways.add((pid, 'B', paid, pname))
                    gids = bcids[gene_id]
                    for gid in gids:
                        gene_pathways.add((gid, pid))

    pathways = list(pathways)
    gene_pathways = list(gene_pathways)

    return pathways, gene_pathways

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config, dbutils

    pathways, gene_pathways = parse(open(sys.argv[1]), open(sys.argv[2]), \
        dbutils.KOBASDB(config.getrc()['kobasdb'] + sys.argv[3] + '.db'))

    print len(pathways), len(gene_pathways)

    g, p = set(), set()
    for gp in gene_pathways:
        g.add(gp[0])
        p.add(gp[1])

    print len(g), len(p)

    pprint(pathways[:5])
    pprint(gene_pathways[:5])
