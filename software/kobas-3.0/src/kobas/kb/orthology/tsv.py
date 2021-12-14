#!/usr/bin/env python
import os

from kobas import dbutils

def get_geid_oeids(orthologydir, abbr):
    geid_oeids = {}

    fns = os.listdir(orthologydir)
#    fns = ['hsa_mmu.tsv']
    for fn in fns:
        species, ospecies = fn.split('.')[0].split('_')
        if species == abbr:
            handle = open(orthologydir + fn)
            for line in handle:
                geid, oeid = line.strip().split('\t')
                geid_oeids.setdefault(geid, set()).add((oeid, ospecies))
        elif ospecies == abbr:
            handle = open(orthologydir + fn)
            for line in handle:
                oeid, geid = line.strip().split('\t')
                geid_oeids.setdefault(geid, set()).add((oeid, species))

    return geid_oeids

def parse(kobasrc, orthologydir, abbr):
    orthologs = set()

    speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + abbr + '.db')
    gids = [gid[0] for gid in speciesdb.genes(name = False)]

    geid_oeids = get_geid_oeids(orthologydir, abbr)

    for gid in gids:
        geids = speciesdb.ensembl_gene_ids_from_gid(gid)
        for geid in geids:
            if geid_oeids.has_key(geid[0]):
                for oeid, ospecies in geid_oeids[geid[0]]:
                    ospeciedb = dbutils.KOBASDB(kobasrc['kobasdb'] + ospecies + '.db')
                    oids = ospeciedb.gids_from_ensembl_gene_id(oeid)
                    for oid in oids:
                        orthologs.add((gid, oid[0]))

    return list(orthologs)

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config

    kobasrc = config.getrc()
    orthologydir = kobasrc['kobas_home'] + '/orthology/'
    orthologs = parse(kobasrc, orthologydir, 'hsa')

    print len(orthologs)
    pprint(orthologs[:5])
