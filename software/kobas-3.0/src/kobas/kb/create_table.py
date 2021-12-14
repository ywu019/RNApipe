#!/usr/bin/env python

def organism(con):
    con.executescript(
        '''
CREATE TABLE Organisms
(
    abbr TEXT PRIMARY KEY,
    name TEXT
);
        ''')

def species(con):
    con.executescript(
        '''
CREATE TABLE Genes
(
    gid TEXT PRIMARY KEY,
    name TEXT
);
CREATE TABLE GeneEntrezGeneIds
(
    gid TEXT,
    entrez_gene_id TEXT,
    PRIMARY KEY (gid, entrez_gene_id)
);
CREATE INDEX GeneEntrezGeneIds_idx_entrez_gene_id ON GeneEntrezGeneIds (entrez_gene_id);
CREATE TABLE GeneGis
(
    gid TEXT,
    gi TEXT,
    PRIMARY KEY (gid, gi)
);
CREATE INDEX GeneGis_idx_gi ON GeneGis (gi);
CREATE TABLE GeneUniprotkbAcs
(
    gid TEXT,
    uniprotkb_ac TEXT,
    PRIMARY KEY (gid, uniprotkb_ac)
);
CREATE INDEX GeneUniprotkbAcs_idx_uniprotkb_ac ON GeneUniprotkbAcs (uniprotkb_ac);
CREATE TABLE GeneEnsemblGeneIds
(
    gid TEXT,
    ensembl_gene_id TEXT,
    PRIMARY KEY (gid, ensembl_gene_id)
);
CREATE INDEX GeneEnsemblGeneIds_idx_ensembl_gene_id ON GeneEnsemblGeneIds (ensembl_gene_id);
CREATE TABLE Orthologs
(
    gid TEXT,
    oid TEXT,
    PRIMARY KEY (gid, oid)
);
CREATE TABLE Pathways
(
    pid INTEGER PRIMARY KEY,
    db TEXT,
    id TEXT,
    name TEXT
);
CREATE TABLE GenePathways
(
    gid TEXT,
    pid INTEGER,
    PRIMARY KEY (gid, pid)
);
CREATE TABLE Diseases
(
    did INTEGER PRIMARY KEY,
    db TEXT,
    id TEXT,
    name TEXT
);
CREATE TABLE GeneDiseases
(
    gid TEXT,
    did INTEGER,
    PRIMARY KEY (gid, did)
);
CREATE TABLE Gos
(
    goid TEXT PRIMARY KEY,
    name TEXT
);
CREATE TABLE GeneGos
(
    gid TEXT,
    goid TEXT,
    PRIMARY KEY (gid, goid)
);
        ''')

def ko(con):
    con.executescript(
        '''
CREATE TABLE Kos
(
    koid TEXT PRIMARY KEY,
    name TEXT
);
CREATE TABLE KoGenes
(
    koid TEXT,
    gid TEXT,
    PRIMARY KEY (koid, gid)
);
CREATE INDEX KoGenes_idx_gid ON KoGenes (gid);
CREATE TABLE KoEntrezGeneIds
(
    koid TEXT,
    entrez_gene_id TEXT,
    PRIMARY KEY (koid, entrez_gene_id)
);
CREATE INDEX KoEntrezGeneIds_idx_entrez_gene_id ON KoEntrezGeneIds (entrez_gene_id);
CREATE TABLE KoGis
(
    koid TEXT,
    gi TEXT,
    PRIMARY KEY (koid, gi)
);
CREATE INDEX KoGis_idx_gi ON KoGis (gi);
CREATE TABLE KoUniprotkbAcs
(
    koid TEXT,
    uniprotkb_ac TEXT,
    PRIMARY KEY (koid, uniprotkb_ac)
);
CREATE INDEX KoUniprotkbAcs_idx_uniprotkb_ac ON KoUniprotkbAcs (uniprotkb_ac);
CREATE TABLE KoEnsemblGeneIds
(
    koid TEXT,
    ensembl_gene_id TEXT,
    PRIMARY KEY (koid, ensembl_gene_id)
);
CREATE INDEX KoEnsemblGeneIds_idx_ensembl_gene_id ON KoEnsemblGeneIds (ensembl_gene_id);
CREATE TABLE Pathways
(
    pid TEXT PRIMARY KEY,
    name TEXT
);
CREATE TABLE KoPathways
(
    koid TEXT,
    pid TEXT,
    PRIMARY KEY (koid, pid)
);
        ''')
