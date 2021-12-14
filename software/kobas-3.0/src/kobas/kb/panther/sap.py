#!/usr/bin/env python
from kobas import dbutils
from kobas.kb.panther import panther_organisms

def parse(handle, kobasrc):
    pathways, gene_pathways = {}, {}

    for line in handle:
        paid, pname, tmp1, tmp2, gene = line.split('\t')[:5]

        sname = gene.split('|')[0]
        if sname not in panther_organisms.panther_organisms.keys():
        # if sname not in ['HUMAN', 'MOUSE']:
            continue
        abbr = panther_organisms.panther_organisms[sname]
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

        if 'UniProtKB' in gene:
            uniprotkb_ac = gene[gene.find('UniProtKB'):].split('|')[0].split('=')[1]
            gids = speciesdb.gids_from_uniprotkb_ac(uniprotkb_ac)
        elif 'ENTREZ' in gene:
            entrez_gene_id = gene[gene.find('ENTREZ'):].split('|')[0].split('=')[1]
            gids = speciesdb.gids_from_entrez_gene_id(entrez_gene_id)
        else:
            continue

        for gid in gids:
            if not pathways.setdefault(abbr, {}).has_key(pname):
                pathways[abbr][pname] = [0, paid]
            gene_pathways.setdefault(abbr, set()).add((gid[0], pname))

    for abbr in pathways.keys():
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')

        count = 0
        for pname in pathways[abbr].keys():
            count += 1
            pathways[abbr][pname][0] = count * 10 + 6

            if __name__ == '__main__':
                print pathways[abbr][pname][0], 'p', pathways[abbr][pname][1], pname
            else:
                speciesdb.con.execute('INSERT OR IGNORE INTO Pathways VALUES (?, ?, ?, ?)', \
                    (pathways[abbr][pname][0], 'p', pathways[abbr][pname][1], pname))

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
