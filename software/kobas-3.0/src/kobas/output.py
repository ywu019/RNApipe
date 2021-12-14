#!/usr/bin/env python

import dbutils, discover

from kobas.kb.biocyc import biocyc_organisms

STOPLIST = set(['of', 'a', 'an', 'to', 'on', 'and', 'with', 'or', 'not', 'syndrome', \
    'by', '&', 'in', 'progression', 'the', 'recurrent', '-', 'disease', 'secondary', \
    'and/or', 'due', 'as', 'syndrome', 'traits', 'chronic', 'familial', 'type', '1', \
    '2', '3', '4', '5', '6', '7', '8', 'primary', 'disorder'])

def annotate_table(items, species):
    for item in items:
        slinks = []
        if len(item.links) > 0:
            for link in item.links:
                if species == 'ko':
                    hyperlink = 'http://www.genome.jp/dbget-bin/www_bget?ko:' + link[0]
                else:
                    hyperlink = 'http://www.genome.jp/dbget-bin/www_bget?' + link[0]
                slink = '|'.join([link[0], link[1], hyperlink])
                slinks.append(slink)
        else:
            slinks.append('None')
        print '%s\t%s' % (item.query, '!'.join(slinks))

def annotate_text(items, species, speciesdb):
    print '\n' + '-' * 20 + '\n'
    print '////'    ##for webserver

    for item in items:
        query = item.query
        print 'Query:' + ' ' * (20 - len('Query:')) + '\t' + query

        if species == 'ko':
            kos = item.links
            if kos:
                pathways = set()
                line_no = 1
                for ko in kos:
                    if line_no == 1:
                        print 'KO:' + ' ' * (20 - len('KO:')) + '\t' + ko[0] + '\t' + ko[1]
                        line_no += 1
                    else:
                        print ' ' * 20 + '\t' + ko[0] + '\t' + ko[1]
                    for pathway in speciesdb.pathways_from_koid(ko[0]):
                        pathways.add('\t'.join([pathway['name'], dbutils.P['K'], pathway['pid']]))
                if pathways:
                    line_no = 1
                    for pathway in pathways:
                        if line_no == 1:
                            print 'Pathway:' + ' ' * (20 - len('Pathway:')) + '\t' + pathway
                            line_no += 1
                        else:
                            print ' ' * 20 + '\t' + pathway
        else:
            if len(item.links) == 1:
                gene = item.links.pop()
                print 'Gene:' + ' ' * (20 - len('Gene:')) + '\t' + gene[0] + '\t' + gene[1]

                entrez_gene_id = speciesdb.entrez_gene_ids_from_gid(gene[0]).fetchone()
                if entrez_gene_id:
                    print 'Entrez Gene ID:' + ' ' * (20 - len('Entrez Gene ID:')) + ' \t' + entrez_gene_id[0]

                pathways = speciesdb.pathways_from_gid(gene[0])
                line_no = 1
                for pathway in pathways:
                    pathway = '\t'.join([pathway['name'], dbutils.P[pathway['db']], pathway['id']])
                    if line_no == 1:
                        print 'Pathway:' + ' ' * (20 - len('Pathway:')) + '\t' + pathway
                        line_no += 1
                    else:
                        print ' ' * 20 + '\t' + pathway

                diseases = speciesdb.diseases_from_gid(gene[0])
                line_no = 1
                for disease in diseases:
                    disease = '\t'.join([disease['name'], dbutils.D[disease['db']], disease['id']])
                    if line_no == 1:
                        print 'Disease:' + ' ' * (20 - len('Disease:')) + '\t' + disease
                        line_no += 1
                    else:
                        print ' ' * 20 + '\t' + disease

                gos = speciesdb.gos_from_gid(gene[0])
                line_no = 1
                for go in gos:
                    go = '\t'.join([go['name'], dbutils.G['G'], go['goid']])
                    if line_no == 1:
                        print 'GO:' + ' ' * (20 - len('GO:')) + '\t' + go
                        line_no += 1
                    else:
                        print ' ' * 20 + '\t' + go

                goslims = speciesdb.goslims_from_gid(gene[0])
                line_no = 1
                for goslim in goslims:
                    goslim = '\t'.join([goslim['name'], dbutils.G['S'], goslim['goid']])
                    if line_no == 1:
                        print 'GOslim:' + ' ' * (20 - len('GOslim:')) + '\t' + goslim
                        line_no += 1
                    else:
                        print ' ' * 20 + '\t' + goslim


        print '////'

def identify(res, abbr, odistr):
    res_p = discover.TestResult(res.title + ['Hyperlink'])
    res_d = discover.TestResult(res.title + ['Hyperlink'])
    res_g = discover.TestResult(res.title + ['Hyperlink'])

    for row in res:
        if row[1] in dbutils.P.keys():
            res_p.add_row([row[0], dbutils.P[row[1]]] + row[2:] + [get_hyperlink(row[1], row[2], abbr, odistr)])
        elif row[1] in dbutils.D.keys():
            res_d.add_row([row[0], dbutils.D[row[1]]] + row[2:] + [get_hyperlink(row[1], row[2])])
        elif row[1] in dbutils.G.keys():
            res_g.add_row([row[0], dbutils.G[row[1]]] + row[2:] + [get_hyperlink(row[1], row[2], abbr)])

    res_p = group(res_p)
    res_d = group(res_d)

    print
    print res_p
    print '-' * 20 + '\n'
    print res_d
    print '-' * 20 + '\n'
    print res_g
    print '-' * 20

def get_hyperlink(dbtype, dbid, abbr = None, odistr = None):
    if dbid:
        if dbtype == 'K':
            url = 'http://www.genome.jp/kegg-bin/show_pathway?' + dbid
            for link in odistr[dbid]:
                url = url + '/' + link + '%09red'
        elif dbtype == 'n':
            url = 'http://pid.nci.nih.gov/search/pathway_landing.shtml?pathway_id=' + dbid + \
                '&source=NCI-Nature%20curated&what=graphic&jpg=on&ppage=1'
        elif dbtype == 'b':
            url = 'http://pid.nci.nih.gov/search/pathway_landing.shtml?pathway_id=' + dbid + \
                '&source=BioCarta%20Imported&what=graphic&jpg=on&ppage=1'
        elif dbtype == 'R':
            url = 'http://www.reactome.org/cgi-bin/eventbrowser_st_id?ST_ID=' + dbid
        elif dbtype == 'B':
            organism = biocyc_organisms.biocyc_organisms2[abbr]
            url = 'http://biocyc.org/%s/NEW-IMAGE?type=NIL&object=%s' % (organism, dbid)
        elif dbtype == 'p':
            url = 'http://www.pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=' + dbid
        elif dbtype == 'o':
            url = 'http://omim.org/entry/' + dbid
        elif dbtype == 'k':
            url = 'http://www.genome.jp/dbget-bin/www_bget?' + dbid
        elif dbtype == 'f':
            url = 'http://django.nubic.northwestern.edu/fundo/databrowse/do_rif/doliteterm/objects/' + dbid
        elif dbtype == 'G':
            url = 'http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=' + dbid
        elif dbtype =='S':
            url = 'http://amigo.geneontology.org/cgi-bin/amigo/term_details?term=' + dbid
    else:
        url = 'None'
    return url

def group(res):
    redun_check, display = [], []
    for row in res:
        name = set([word for word in row[0].lower().split() if word not in STOPLIST])
        if name not in redun_check:
            redun_check.append(name)
            display.append(row)
        else:
            index = redun_check.index(name)
            redun_check.insert(index, set())
            display.insert(index + 1, ['@@' + row[0]] + row[1:])
    result = discover.TestResult(res.title)
    for row in display:
        result.add_row(row)
    return result
