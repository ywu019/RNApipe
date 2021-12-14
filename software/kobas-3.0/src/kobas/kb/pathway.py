#!/usr/bin/env python

from kobas.kb.pid import pid_xml
from kobas.kb.reactome import reactome
from kobas.kb.panther import sap
from kobas.kb.biocyc import col

def insert_pid(kobasdb, kobasdir):
    ppathways, pgene_pathways = pid_xml.parse(open(kobasdir + '/pid/NCI-Nature_Curated.xml'), kobasdb)
    kobasdb.con.executemany('INSERT OR IGNORE INTO Pathways VALUES (?, ?, ?, ?)', ppathways)
    kobasdb.con.executemany('INSERT OR IGNORE INTO GenePathways VALUES (?, ?)', pgene_pathways)
    bpathways, bgene_pathways = pid_xml.parse(open(kobasdir + '/biocarta/BioCarta.xml'), kobasdb)
    kobasdb.con.executemany('INSERT OR IGNORE INTO Pathways VALUES (?, ?, ?, ?)', bpathways)
    kobasdb.con.executemany('INSERT OR IGNORE INTO GenePathways VALUES (?, ?)', bgene_pathways)

def insert_reactome(kobasrc, kobasdir):
    reactome.parse(open(kobasdir + '/reactome/UniProt2Reactome_All_Levels.txt'), \
        kobasrc)

def insert_panther(kobasrc, kobasdir):
    sap.parse(open(kobasdir + '/panther/SequenceAssociationPathway3.4.txt'), kobasrc)

def insert_biocyc(kobasdb, kobasdir, bcdir):
    pathways, gene_pathways = col.parse(open(kobasdir + '/biocyc/' + bcdir + '/genes.col'), \
            open(kobasdir + '/biocyc/' + bcdir + '/pathways.col'), kobasdb)
    kobasdb.con.executemany('INSERT OR IGNORE INTO Pathways VALUES (?, ?, ?, ?)', pathways)
    kobasdb.con.executemany('INSERT OR IGNORE INTO GenePathways VALUES (?, ?)', gene_pathways)
