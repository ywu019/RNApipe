#!/usr/bin/env python

from kobas import dbutils
from kobas.kb.go import obo, go_organisms, association

go_categories = ('GO:0003674', 'GO:0005575', 'GO:0008150')
##molecular_function, cellular_component, biological_process

def insert_go(kobasrc, kobasdir):
    go_terms, go_inferreds = obo.parse(open(kobasdir + '/go/gene_ontology.1_2.obo'))

    for organism in go_organisms.go_organisms:
#    for organism in (('gene_association.sgd', 'gp2protein.sgd', '559292', 'sce'), \
#          ('gene_association.ecocyc', 'not_need', '83333', 'eco')):
        print organism
    #for organism in (('gene_association.goa_human', 'not_need', '9606', 'hsa'), \
    #     ('gene_association.mgi', 'gp2protein.mgi', '10090', 'mmu')):
        multi_uniprotkb_ac_gos = association.parse(kobasdir, organism)

        abbrs = [abbr for abbr in organism[3].split('|')]
        for abbr in abbrs:
            speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')
            gos, gene_gos = set(), set()
            for uniprotkb_ac_go in multi_uniprotkb_ac_gos[abbr]:
                if (uniprotkb_ac_go[1] not in go_categories) and (uniprotkb_ac_go[1] in go_terms.keys()):
                    gids = speciesdb.gids_from_uniprotkb_ac(uniprotkb_ac_go[0])
                    for gid in gids:
                        gene_gos.add((gid[0], uniprotkb_ac_go[1]))
                        gos.add((uniprotkb_ac_go[1], go_terms[uniprotkb_ac_go[1]]))
                        for goid in go_inferreds[uniprotkb_ac_go[1]]:
                            if (goid not in go_categories) and (goid in go_terms.keys()):
                                gene_gos.add((gid[0], goid))
                                gos.add((goid, go_terms[goid]))

            speciesdb.con.executemany('INSERT OR IGNORE INTO Gos VALUES (?, ?)', list(gos))
            speciesdb.con.executemany('INSERT OR IGNORE INTO GeneGos VALUES (?, ?)', list(gene_gos))

