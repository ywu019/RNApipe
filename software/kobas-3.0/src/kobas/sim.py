#!/usr/bin/env python
import os, time, math
import numpy as np
#import rpy2.robjects as robjects
from itertools import combinations
from kobas import discover, annot, dbutils, config 

FILENAME = """/rd1/user/tangkj/trunk/test/result_I"""
kobasrc = config.getrc()
speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + 'hsa' + '.db')

COLUMN = ('Term', 'Database', 'ID', 'Input number', 'Background number', 'P-Value', 'Corrected P-Value', 'Input', 'Hyperlink')
#robjects.r.library('epiR')

def verify_file(file_name):
    ##verify file existence and accession, and return file handle
    if os.access(file_name, os.F_OK):
        return open(file_name)
    else:
        print 'file %s does not exist' % file_name
        sys.exit(1)

def oneTerm(line):
    cont = line.split('\t')
    term = dict(zip(COLUMN, cont))
    return term
    
def cut(handle):
    Plist =[]
    Dlist = []
    Glist = []    
    flag = 0
    for line in handle.readlines():
        if line.startswith('#Term'):
            flag += 1
        if line.strip() and not line.startswith('#Term') and not line.startswith('-'):
            if flag == 1:
                Plist.append(oneTerm(line.strip()))
            if flag == 2:
                Dlist.append(oneTerm(line.strip()))
            if flag == 3:
                Glist.append(oneTerm(line.strip()))
    return Plist, Dlist, Glist
"""
def makeMatrix(termA, termB):
    matrix = [0, 0, 0, 0]
    for i in range(len(termA)):
        if termA[i] == termB[i] == 0:
            matrix[3] += 1
        elif termA[i] == termB[i] == 1:
            matrix[0] += 1
        elif termA[i] == 1 and termB[i] == 0:
            matrix[1] += 1
        else:
            matrix[2] += 1
    return matrix
"""
def makeMatrix(termA, termB, num):
    matrix = [0, 0, 0, 0]
    lenA = len(termA)
    lenB = len(termB)
    '''
    i, j = 0, 0
    while(i < lenA and j < lenB):
	if termA[i] < termB[j]:
	    matrix[1] += 1
	    i += 1
	elif termA[i] == termB[j]:
            matrix[0] += 1
            i += 1
	    j += 1
	else:
	    matrix[2] += 1
            j += 1
    if i == lenA:
	matrix[2] += (lenB - j)
    if j == lenB:
	matrix[1] += (lenA - i)
    matrix[3] = num - sum(matrix[:3])
    return matrix'''
    matrix[0] = len(set(termA) & set(termB))
    matrix[1] = lenA - matrix[0]
    matrix[2] = lenB - matrix[0]
    matrix[3] = num - sum(matrix[:3])
    return matrix
         
def kappa(termA, termB, num):
    matrix = makeMatrix(termA, termB, num)
    #Tab = sum(matrix)
    Tab = num
    Oab = (matrix[0] + matrix[3]) / (Tab + 0.0)
    Aab = ((matrix[0] + matrix[1]) * (matrix[0] + matrix[2]) + (matrix[2] + matrix[3]) * (matrix[1] + matrix[3])) / (Tab * Tab + 0.0)
    Kab = (Oab - Aab) / (1.0 - Aab)
    #cmd = """epi.kappa(matrix(c(%d, %d, %d, %d), nc = 2))"""%tuple(matrix)
    #Kab = robjects.r(cmd)[0][0][0]
    return Kab

def Jaccard(termA, termB, num):
    numerator = 0.0
    denominator = 0.0
    '''
    for i in range(len(termA)):
        if termA[i] == termB[i] == 1:
            numerator += 1
            denominator += 1
        elif termA[i] == termB[i] == 0:
            pass
        else:
            denominator += 1'''
    numerator += len(set(termA) & set(termB))
    denominator += len(set(termA) | set(termB))
    return numerator / denominator

def termNum(Plist, Dlist, Glist):
    vector = []
    for term in (Plist + Dlist + Glist):
        vector.extend(term['Input'].split('|'))
    return len(set(vector))

"""
def initVector(Plist, Dlist, Glist):
    vector = []
    for term in (Plist + Dlist + Glist):
        vector.extend(term['Input'].split('|'))
    return list(set(vector))


def termVector(term, vector):
    termset = term['Input'].split('|')
    termvector = []
    for i in vector:
        if i in termset:
            termvector.append(1)
        else:
            termvector.append(0)
    return termvector
"""

def termVector(term):
    termset = sorted(term['Input'].split('|'))
    return termset
"""
def listMatrix(Xlist, vector, method):
    matrix = []
    for i in range(len(Xlist)):
        matrix.append([])
        vectorA = termVector(Xlist[i], vector)
        for j in range(len(Xlist)):
            vectorB = termVector(Xlist[j], vector)
            matrix[i].append(method(vectorA, vectorB))
    return matrix
"""
def listMatrix(Xlist, num, method):
    size = len(Xlist)
    matrix = np.empty((size, size), np.float)
    for i in range(size):
        termA = termVector(Xlist[i])
        for j in range(i+1, size):
            termB = termVector(Xlist[j])
            matrix[i][j] = method(termA, termB, num)
            matrix[j][i] = matrix[i][j]
    return list(matrix)

    
def testSelf(matrix, clus, thresh):
    comb = combinations(clus, 2)
    count = 0.0
    ok = 0
    for i in comb:
        count += 1
        if matrix[i[0]][i[1]] >= thresh:
            ok += 1
    if count > 3 and (ok / count) > 0.5:
        return True
    else:
        return False

def testTwo(clus1, clus2):
    len1 = len(clus1)
    len2 = len(clus2)
    #count = 0.0
    two = list(set(clus1) | set(clus2))
    
    if len(two) < (0.5 * max(len1, len2) + min(len1, len2)):
        return two
    else:
        return False
    '''
    for i in clus1:
        if i in clus2:
            count += 1
    if (count / len1) > 0.5 and (count / len2) > 0.5:
        return two
    else:
        return False'''

def clusterInit(matrix, thresh):
    allclus = []
    for i in range(len(matrix)):
        clus = [i]
        for j in range(len(matrix)):
            if j != i and matrix[i][j] >= thresh:
                clus.append(j)
        if testSelf(matrix, clus, thresh):
            allclus.append(clus)
    return allclus
"""
def clusterIter(matrix, allclus):
    #print allclus
    for i in range(len(allclus)):
        for j in range(i+1, len(allclus)):
            clusA, clusB = allclus[i], allclus[j]
            two =  testTwo(clusA, clusB)
                #allclus.append(list(set(clusA) | set(clusB)))
                #allclus = allclus[:i] + allclus[i+1:j] + allclus[j+1:] + [list(set(clusA) | set(clusB))]
            if two:     
                allclus = [two] + allclus 
                allclus.remove(clusA)
                allclus.remove(clusB)
                return clusterIter(matrix, allclus)
    return allclus
"""
def clusterIter(matrix, allclus):
    i = 0
    j = 1
    lenL = len(allclus)
    while(i < lenL - 1):
        clusA, clusB = allclus[i], allclus[j]
        two = testTwo(clusA, clusB)
        if two:
            allclus.append(two)
            allclus.remove(clusA)
            allclus.remove(clusB)
            lenL -= 1
            j = i + 1
        else:
            if j < lenL - 1:
                j += 1
            else:
                i += 1
                j = i + 1
    return allclus

def calculateScore(clus, Xlist):
    mul = 1.0
    for i in clus:
        mul *= (-math.log10(float(Xlist[i]['P-Value'])))
    return mul ** (1./len(clus))
        
if __name__ == "__main__":
    '''
    handle = verify_file(FILENAME)
    P, D, G = cut(handle)
    vector = initVector(P, D, G)
    begin = time.time()
    matrix = listMatrix(P, vector, kappa)
    print len(matrix)
    print time.time() - begin
    '''
    matrix = [[1.0,1.0,1.0,0.35,-0.50,-0.50,-0.50,0.00],
              [1.0,1.0,1.0,0.35,-0.50,-0.50,-0.50,0.00],
              [1.0,1.0,1.0,0.35,-0.50,-0.50,-0.50,0.00],
              [0.35,0.35,0.35,1.0,0.35,0.35,0.35,-0.11],
              [-0.50,-0.50,-0.50,0.35,1.0,1.0,1.0,0.00],
              [-0.50,-0.50,-0.50,0.35,1.0,1.0,1.0,0.00],
              [-0.50,-0.50,-0.50,0.35,1.0,1.0,1.0,0.00],
              [0.00,0.00,0.00,-0.11,0.00,0.00,0.00,1.0]]
    
    #begin = time.time()
    allclus = clusterInit(matrix, 0.35)
    allclus = clusterIter(matrix, allclus)
    '''
    print time.time() - begin
    print len(allclus)
    #print allclus
    for i in allclus[0]:
        print P[i]'''
