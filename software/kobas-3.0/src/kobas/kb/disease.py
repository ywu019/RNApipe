#!/usr/bin/env python

from kobas.kb.omim import morbidmap
from kobas.kb.fundo import do_lite
from kobas.kb.gad import gad_all
from kobas.kb.nhgri import gwascatalog

def insert_omim(kobasdb, kobasdir):
    diseases, gene_diseases = morbidmap.parse(open(kobasdir + '/omim/morbidmap'), \
        open(kobasdir + '/omim/human_gene_info'), kobasdb)
    kobasdb.con.executemany('INSERT INTO Diseases VALUES (?, ?, ?, ?)', diseases)
    kobasdb.con.executemany('INSERT INTO GeneDiseases VALUES (?, ?)', gene_diseases)

def insert_fundo(kobasdb, kobasdir):
    diseases, gene_diseases = do_lite.parse(open(kobasdir + '/fundo/do_lite.txt'), \
        open(kobasdir + '/fundo/DO_Lite_Terms.html'), kobasdb)
    kobasdb.con.executemany('INSERT INTO Diseases VALUES (?, ?, ?, ?)', diseases)
    kobasdb.con.executemany('INSERT INTO GeneDiseases VALUES (?, ?)', gene_diseases)

def insert_gad(kobasdb, kobasdir):
    diseases, gene_diseases = gad_all.parse(open(kobasdir + '/gad/all.txt'), kobasdb)
    kobasdb.con.executemany('INSERT INTO Diseases VALUES (?, ?, ?, ?)', diseases)
    kobasdb.con.executemany('INSERT INTO GeneDiseases VALUES (?, ?)', gene_diseases)

def insert_nhgri(kobasdb, kobasdir):
    diseases, gene_diseases = gwascatalog.parse(open(kobasdir + '/nhgri/gwascatalog'), \
        open(kobasdir + '/omim/human_gene_info'), kobasdb)
    kobasdb.con.executemany('INSERT INTO Diseases VALUES (?, ?, ?, ?)', diseases)
    kobasdb.con.executemany('INSERT INTO GeneDiseases VALUES (?, ?)', gene_diseases)
