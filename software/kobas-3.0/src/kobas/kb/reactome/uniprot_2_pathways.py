#!/usr/bin/env python

import re

from kobas import dbutils
from kobas.kb.reactome import reactome_organisms

def parse(stid_handle, map_handle, kobasrc):
    pathways, gene_pathways = {}, {}

    stids = {}
    for line in stid_handle:
        uniprotkb_ac, paid, pname = line.split('\t')[:3]
        if not stids.has_key((uniprotkb_ac, pname)):
            stids[(uniprotkb_ac, pname)] = paid

    for line in map_handle:
        uniprotkb_ac, tmp1, pnames, tmp2, sname = line[:-1].split('\t')

        if sname not in reactome_organisms.reactome_organisms.keys():
        # if sname not in ['Homo sapiens', 'Mus musculus']:
            continue
        abbr = reactome_organisms.reactome_organisms[sname]
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

        gids = speciesdb.gids_from_uniprotkb_ac(uniprotkb_ac)
        pnames = re.sub('\[.*\]: ', '', pnames).split('; ')
        for gid in gids:
            for pname in pnames:
                pname = re.sub(' IEA$', '', pname)
                if stids.has_key((uniprotkb_ac, pname)):
                    paid = stids[(uniprotkb_ac, pname)]
                else:
                    paid = ''

                if not pathways.setdefault(abbr, {}).has_key(pname):
                    pathways[abbr][pname] = [0, paid]
                gene_pathways.setdefault(abbr, set()).add((gid[0], pname))

    for abbr in pathways.keys():
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

        count = 0
        for pname in pathways[abbr].keys():
            count += 1
            pathways[abbr][pname][0] = count * 10 + 4

            if __name__ == '__main__':
                print pathways[abbr][pname][0], 'R', pathways[abbr][pname][1], pname
            else:
                speciesdb.con.execute('INSERT INTO Pathways VALUES (?, ?, ?, ?)', \
                    (pathways[abbr][pname][0], 'R', pathways[abbr][pname][1], pname))

        for gene_pathway in gene_pathways[abbr]:
            if __name__ == '__main__':
                print gene_pathway[0], pathways[abbr][gene_pathway[1]][0]
            else:
                speciesdb.con.execute('INSERT INTO GenePathways VALUES (?, ?)', \
                    (gene_pathway[0], pathways[abbr][gene_pathway[1]][0]))

if __name__ == '__main__':
    import sys

    from kobas import config

    parse(open(sys.argv[1]), open(sys.argv[2]), config.getrc())
