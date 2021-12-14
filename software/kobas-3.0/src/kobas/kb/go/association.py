#!/usr/bin/env python

def parse(kobasdir, organism):
    multi_uniprotkb_ac_gos = {}

    if organism[1] != 'not_need':
        gp2protein_handle = open(kobasdir + '/go/gp2protein/' + organism[1])
        gp2protein = {}
        for line in gp2protein_handle:
            if line[0] != '!':
                info = line[:-1].split('\t')
                if len(info) == 2:
                    gpid = info[0].split(':')[-1]
                    index = info[1].find('UniProtKB:')
                    if index != -1:
                        uniprotkb_ac = info[1][index + 10:index + 16]
                        gp2protein.setdefault(gpid, set()).add(uniprotkb_ac)

    tids = [tid for tid in organism[2].split('|')]
    abbrs = [abbr for abbr in organism[3].split('|')]
    tid_abbrs = {}
    for index in range(len(tids)):
        tid_abbrs[tids[index]] = abbrs[index]
        multi_uniprotkb_ac_gos[abbrs[index]] = set()

    association_handle = open(kobasdir + '/go/association/' + organism[0])
    association = {}
    for line in association_handle:
        if line[0] != '!':
            info = line[:-1].split('\t')
            if 'not' not in info[3].lower():    ##if 'not' in column 4, discard the annotation
                tid = info[12].split('|')[0][6:]    ##if '|' in column 13, only the first organism encodes the gene or gene product
                if tid_abbrs.has_key(tid):
                    abbr = tid_abbrs[tid]
                    gpid = info[1].split(':')[-1]    ##column 2 is DB Object ID
                    if organism[1] != 'not_need':
                        uniprotkb_acs = gp2protein.get(gpid, [])
                    else:
                        uniprotkb_acs = [gpid]
                    goid = info[4]    ##column 5 is GO ID
                    for uniprotkb_ac in uniprotkb_acs:
                        multi_uniprotkb_ac_gos[abbr].add((uniprotkb_ac, goid))

    return multi_uniprotkb_ac_gos

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config

    multi_uniprotkb_ac_gos = parse(config.getrc()['kobas_home'], ('gene_association.ecocyc', 'not_need', '83333', 'eco'))
#    multi_uniprotkb_ac_gos = parse(config.getrc()['kobas_home'], ('gene_association.goa_human', 'not_need', '9606', 'hsa'))

    print len(multi_uniprotkb_ac_gos)
    print len(multi_uniprotkb_ac_gos['eco'])
    pprint(tuple(multi_uniprotkb_ac_gos['eco'])[:5])
