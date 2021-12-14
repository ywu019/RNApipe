#!/usr/bin/env python

from kobas.kb.kegg import ko, ent, br08402_gene

def insert_organism(con, kobasdir):
    organism_file = open(kobasdir + '/kegg/avai_taxonomy')
    for line in organism_file:
        abbr, tmp1, name = line[:-1].split('\t')
        con.execute('INSERT INTO Organisms VALUES (?, ?)', (abbr, name))

def insert_ko(con, kobasdir):
    kos, pathways, ko_pathways = ko.parse(open(kobasdir + '/kegg/ko'))
    con.executemany('INSERT INTO Kos VALUES (?, ?)', kos)
    con.executemany('INSERT INTO Pathways VALUES (?, ?)', pathways)
    con.executemany('INSERT INTO KoPathways VALUES (?, ?)', ko_pathways)

def insert_ent(spcon, kocon, kobasdir, abbr):
    genes, pathways, gene_pathways, ko_genes, \
    gene_entrez_gene_ids, gene_gis, gene_uniprotkb_acs, gene_ensembl_gene_ids, \
    ko_entrez_gene_ids, ko_gis, ko_uniprotkb_acs, ko_ensembl_gene_ids \
    = ent.parse(open(kobasdir + '/kegg/organisms/' + abbr + '.ent'), abbr)

    spcon.executemany('INSERT INTO Genes VALUES (?, ?)', genes)
    spcon.executemany('INSERT INTO Pathways VALUES (?, ?, ?, ?)', pathways)
    spcon.executemany('INSERT INTO GenePathways VALUES (?, ?)', gene_pathways)
    spcon.executemany('INSERT INTO GeneEntrezGeneIds VALUES (?, ?)', gene_entrez_gene_ids)
    spcon.executemany('INSERT INTO GeneGis VALUES (?, ?)', gene_gis)
    spcon.executemany('INSERT INTO GeneUniprotkbAcs VALUES (?, ?)', gene_uniprotkb_acs)
    spcon.executemany('INSERT INTO GeneEnsemblGeneIds VALUES (?, ?)', gene_ensembl_gene_ids)
    kocon.executemany('INSERT INTO KoGenes VALUES (?, ?)', ko_genes)
    kocon.executemany('INSERT INTO KoEntrezGeneIds VALUES (?, ?)', ko_entrez_gene_ids)
    kocon.executemany('INSERT INTO KoGis VALUES (?, ?)', ko_gis)
    kocon.executemany('INSERT OR IGNORE INTO KoUniprotkbAcs VALUES (?, ?)', ko_uniprotkb_acs)  ##genes from different species have same UniProtKB ACs in KEGG, which is ugly
    kocon.executemany('INSERT INTO KoEnsemblGeneIds VALUES (?, ?)', ko_ensembl_gene_ids)

def insert_disease(con, kobasdir):
    diseases, gene_diseases = br08402_gene.parse(open(kobasdir + '/kegg/br08402_gene.keg'))
    con.executemany('INSERT INTO Diseases VALUES (?, ?, ?, ?)', diseases)
    con.executemany('INSERT INTO GeneDiseases VALUES (?, ?)', gene_diseases)
