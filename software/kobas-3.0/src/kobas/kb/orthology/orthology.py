#!/usr/bin/env python
from kobas import dbutils
from kobas.kb.orthology import tsv

def insert_orthology(kobasrc, kobasdir):
    organismdb = dbutils.KOBASDB(kobasrc['kobasdb'] + 'organism.db')
    for abbr in organismdb.organisms(name = False):
#    for abbr in [('tru', ), ('xma', )]:
        orthologs = tsv.parse(kobasrc, kobasdir + '/orthology/', abbr[0])
        speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr[0] + '.db')
        speciesdb.con.executemany('INSERT INTO Orthologs VALUES (?, ?)', orthologs)
