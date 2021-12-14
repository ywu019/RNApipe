#!/usr/bin/env python

import re

from kobas import dbutils
from kobas.kb.reactome import reactome_organisms

def parse(handle, kobasrc):
    pathways, gene_pathways = {}, {}

    for line in handle:
        uniprotkb_ac, paid, url, pname, ec, sname = line[:-1].split('\t')
        if sname not in reactome_organisms.reactome_organisms.keys():
        #if sname not in ['Homo sapiens', 'Mus musculus']:
            continue
        abbr = reactome_organisms.reactome_organisms[sname]
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

        gids = speciesdb.gids_from_uniprotkb_ac(uniprotkb_ac)
        for gid in gids:
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
                speciesdb.con.execute('INSERT OR IGNORE INTO Pathways VALUES (?, ?, ?, ?)', \
                    (pathways[abbr][pname][0], 'R', pathways[abbr][pname][1], pname))

        for gene_pathway in gene_pathways[abbr]:
            if __name__ == '__main__':
                print gene_pathway[0], pathways[abbr][gene_pathway[1]][0]
            else:
                speciesdb.con.execute('INSERT OR IGNORE INTO GenePathways VALUES (?, ?)', \
                    (gene_pathway[0], pathways[abbr][gene_pathway[1]][0]))

if __name__ == '__main__':
    import sys

    from kobas import config

    parse(open(sys.argv[1]), config.getrc())
