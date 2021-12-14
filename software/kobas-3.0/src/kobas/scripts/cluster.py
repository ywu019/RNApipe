#!/usr/bin/python

import os, sys, time
from optparse import OptionParser

from kobas import sim

COLUMN = ('Term', 'Database', 'ID', 'Input number', 'Background number', 'P-Value', 'Corrected P-Value', 'Input', 'Hyperlink')

def config_option():
    usage = 'Usage: %prog -i infile [-o outfile] [-t threshold] [-d database] [-m method]'
    p = OptionParser(usage)
    p.add_option(
        '-i', '--infile', dest = 'infile', action = 'store',
        type = 'string', help = 'input file, the output of identify')
    p.add_option(
        '-m', '--method', dest = 'method', default = 'k', action = 'store',
        type = 'string', help = 'clustering method, k for Cohen\'s kappa coefficient, j for Jaccard similarity coefficient')
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for identification result, default stdout')
    p.add_option(
        '-d', '--database', dest = 'data', default = 'P/D/G', action = 'store',
        type = 'string', help = 'databases for selection, 1-letter abbreviation separated by "/": P for PATHWAY, D for DISEASE, G for GO terms, default P/D/G')
    p.add_option(
        '-t', '--threshold', dest = 'threshold', default = 0.35, action = 'store',
        type = 'float', help = 'kappa threshold of clustering, default 0.35')
    opt, args = p.parse_args()
    return (p, opt, args)

opt_parser, opt, args = config_option()

def verify_file(file_name):
    ##verify file existence and accession, and return file handle
    if os.access(file_name, os.F_OK):
        return open(file_name)
    else:
        print 'file %s does not exist' % file_name
        sys.exit(1)

def print_output(X, Xmatrix, threshold):
#    print '-'*20
#    initTime = time.time()
    allclus = sim.clusterInit(Xmatrix, threshold)
#    print 'initTime: %s'%(time.time() - initTime)
#    iterTime = time.time()
    allclus = sim.clusterIter(Xmatrix, allclus)
#    print 'iterTime: %s'%(time.time()-iterTime)
#    sortTime = time.time()
    allclus.sort(key = lambda x:sim.calculateScore(x, X), reverse = True)
#    print 'sortTime: %s'%(time.time()-sortTime)
    for i in range(len(allclus)):
        print '##Cluster: %d\tScore: %s'%(i + 1, sim.calculateScore(allclus[i], X))
        print '#' + '\t'.join(COLUMN)
        allclus[i].sort(key = lambda x:float(X[x]['P-Value']))
        for term in allclus[i]:
            para = []
            for j in COLUMN:
                para.append(X[term][j])
            print '\t'.join(para)
#    print '-'*20

##process opt.fgfile
if not opt.infile:
    opt_parser.error('Option -f must be assigned.\n')
else:
    fg_handle = verify_file(opt.infile)

##process opt.outfile
if opt.outfile:
    global old_stdout
    old_stdout = sys.stdout
    sys.stdout = open(opt.outfile, 'w')

P, D, G = sim.cut(fg_handle)
num = sim.termNum(P, D, G)
dbs=opt.data.split('/')
def cluster(i): 
#    matrixTime = time.time()
    if opt.method == 'k':
        matrix = sim.listMatrix(i, num, sim.kappa)
    if opt.method == 'j':
        matrix = sim.listMatrix(i, num, sim.Jaccard)
#    print 'matrixTime: %s'%(time.time()-matrixTime)
    print_output(i, matrix, opt.threshold)

if 'P' in dbs:
    cluster(P)
if 'D' in dbs:
    cluster(D)
if 'G' in dbs:
    cluster(G)
