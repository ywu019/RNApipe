#!/usr/bin/env python

import sys
#from optparse import OptionParser
from kobas import annot, dbutils, output, config
import pdb

P = {'K': 'KEGG PATHWAY', 'n': 'PID', 'b': 'BioCarta', 'R': 'Reactome', 'B': 'BioCyc', 'p': 'PANTHER'}
#     1                    2           3                4                5              6
D = {'o': 'OMIM', 'k': 'KEGG DISEASE', 'f': 'FunDO', 'g': 'GAD', 'N': 'NHGRI GWAS Catalog'}
#     1            2                    3             4           5
G = {'G': 'Gene Ontology', 'S':'Gene Ontology Slim', 's':'Gene Ontology Slim leaf term'}

DBLINKS = {'ncbigene':'entrez_gene_id', 'ncbigi': 'gi', 'uniprot': 'uniprotkb_ac', 'ensembl': 'ensembl_gene_id'}

kobasrc = config.getrc()

def getGmt(species,dbtype, db, idtype, outdir):
    speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + species + '.db')
    av_dbs = speciesdb.databases_from_abbr(species)
    av_dbs = [each[0] for each in av_dbs]
    if not db in av_dbs :
        print 'Available databases for %s:' % species
        for x in av_dbs:
            print '\t'.join(x)
        sys.exit()
    outf = open(outdir+'/' + ('.').join([species,db,idtype,'gmt']), 'w')
    outf.write('# Species: ' + species + '\n')
    if dbtype == 'P':
        gsets = speciesdb.allpathways(db).fetchall()
#        outf.write('-'*20+'\n')
        outf.write('# Database: '  + P[db] +'\n' )
        for each in gsets:
            genes = []
            pid, Gid, Gname = each
            kgenes = speciesdb.genes_from_pid(pid)
            for kgene in kgenes:
                #speciesdb.dblink_ids_from_gid(kgene[0], idtype).fetchone[0]
                dblink_ids = [ dblink_id[0] for dblink_id in speciesdb.dblink_ids_from_gid(kgene[0], idtype)]
                if(len(dblink_ids)!=0):
                    genes.append(dblink_ids[0])
            outf.write(('\t').join([Gid,Gname]+genes)+'\n')

    if dbtype =='D':
        gsets = speciesdb.alldiseases(db).fetchall()
#        outf.write('-'*20+'\n')
        outf.write('# Database: '  + D[db] +'\n')
        for each in gsets:
            genes = []
            did, Gid, Gname = each
            kgenes = speciesdb.genes_from_did(did)
            for kgene in kgenes:
                dblink_ids = [ dblink_id[0] for dblink_id in speciesdb.dblink_ids_from_gid(kgene[0], idtype)]
                if(len(dblink_ids)!=0):
                    genes.append(dblink_ids[0])
            outf.write(('\t').join([Gid,Gname]+genes)+'\n')
    
    if dbtype == 'G':
        if db == 'G':
            gsets = speciesdb.allgoterms().fetchall()
        if db == 'S':
            gsets = speciesdb.allgoslimterms().fetchall()
#        outf.write('-'*20+'\n')
        outf.write('# Database: ' + G[db] +'\n')
        for each in gsets:
            genes = []
            goid, Gname = each
            kgenes = speciesdb.genes_from_goid(goid)
            for kgene in kgenes:
                dblink_ids = [ dblink_id[0] for dblink_id in speciesdb.dblink_ids_from_gid(kgene[0], idtype)]
                if(len(dblink_ids)!=0):
                    genes.append(dblink_ids[0])
            outf.write(('\t').join([goid,Gname]+genes)+'\n')

    outf.close()





if __name__ == "__main__":
    if(len(sys.argv) != 5) :
        print "Usage: python "+ sys.argv[0] +' species[eg: hsa] db-key-value[eg: P:K] idtype[eg: ensembl/ncbigene/ncbigi/uniprot] output_dir \n'
        sys.exit()
    species = sys.argv[1]
    dbkv = sys.argv[2]
    outdir = sys.argv[4]
    idin = sys.argv[3]
    idtype = DBLINKS[idin]
    dbtype,db = dbkv.split(':')

    getGmt(species,dbtype,db,idtype,outdir)

