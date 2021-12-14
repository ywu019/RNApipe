import sys,os
from kobas import config, dbutils, exception, annot

import rpy2.robjects as robjects



class Distr(dict):

    def __init__(self):
        dict.__init__(self)

    def __getitem__(self, key):
        return dict.get(self, key, set())

    def __setitem__(self, key, val):
        if self.has_key(key):
            self[key].add(val)
        else:
            dict.__setitem__(self, key, set([val]))
#    def size(self,gnum):
#        return gnum

#    def size(self):
#        return sum([len(vals) for vals in self.values()])
#    def size(self):
#        x = set()
#        for vals in self.values():
#            x = x | vals
#            print x 
#        return len(x)

class DistrGnum(object):  
    def __init__(self,distr,gnum = None, speciesdb = None, idmapping = None):
        self.dis = distr
        if gnum:
            self.gnum = int(gnum)
        if (not gnum and speciesdb and idmapping):
            self.gnum = int(speciesdb.gnum_from_dbtype(idmapping))

def remove_small_terms(distr, small_num):
    for term in distr.keys():
        if len(distr[term]) < small_num:
            distr.pop(term)

def remove_small_terms2(distr, small_num):
    for term in distr.keys():
        if distr[term] < small_num:
            distr.pop(term)


##two functions for getting distr
def distr_from_file(items, db, speciesdb, abbr, small_num = None):
    distr, odistr = Distr(), Distr()
    if abbr == 'ko':
        for item in items:    ##item[0] is query, item[1] is links
            for koid in item[1]:
                terms = speciesdb.pathways_from_koid(koid)
                for term in terms:
                    distr[(term['name'], db, term['pid'])] = item[0]
                    odistr[term['pid']] = koid

    elif db in dbutils.P.keys():
        for item in items:
            terms = speciesdb.pathways_from_gid(item[1][0], db)
            for term in terms:
                distr[(term['name'], db, term['id'])] = item[0]
                if db == 'K':
                    odistr[term['id']] = item[1][0]

    elif db in dbutils.D.keys():
        for item in items:
            terms = speciesdb.diseases_from_gid(item[1][0], db)
            for term in terms:
                distr[(term['name'], db, term['id'])] = item[0]

    elif db in dbutils.G.keys():
        for item in items:
            terms = speciesdb.gos_from_gid(item[1][0])
            for term in terms:
                distr[(term['name'], db, term['goid'])] = item[0]

    if small_num:
        remove_small_terms(distr, small_num)
        return distr
    elif db == 'K':
        return distr, odistr
    else:
        return distr

def distr_from_default(default, db, speciesdb, abbr, small_num, idmapping):
    distr = {}
    if abbr == 'ko':
        koids = speciesdb.kos()
        for koid in koids:
            terms = speciesdb.pathways_from_koid(koid['koid'])
            gids = speciesdb.gids_from_koid(koid['koid'], default)
            for gid in gids:
                if idmapping:
                    bspeciesdb = dbutils.KOBASDB(config.getrc()['kobasdb'] + default + '.db')
                    dblink_ids = [dblink_id[0] for dblink_id in bspeciesdb.dblink_ids_from_gid(gid['gid'], idmapping)]
                for term in terms:
                    if idmapping:
                        distr[(term['name'], db, term['pid'])] =  distr.get((term['name'], db, term['pid']), 0) + len(dblink_ids)
                    else:
                        distr[(term['name'], db, term['pid'])] = distr.get((term['name'], db, term['pid']), 0) + 1

    else:
        gids = speciesdb.genes()
        for gid in gids:
            if db in dbutils.P.keys():
                terms = speciesdb.pathways_from_gid(gid['gid'], db)
            elif db in dbutils.D.keys():
                terms = speciesdb.diseases_from_gid(gid['gid'], db)
            elif db in dbutils.G.keys():
                terms = speciesdb.gos_from_gid(gid['gid'])

            if idmapping:
                dblink_ids = [dblink_id[0] for dblink_id in speciesdb.dblink_ids_from_gid(gid['gid'], idmapping)]
            if db in dbutils.P.keys() + dbutils.D.keys():
                for term in terms:
                    if idmapping:
                        distr[(term['name'], db, term['id'])] = distr.get((term['name'], db, term['id']), 0) + len(dblink_ids)
                    else:
                        distr[(term['name'], db, term['id'])] = distr.get((term['name'], db, term['id']), 0) + 1

            elif db in dbutils.G.keys():
                for term in terms:
                    if idmapping:
                        distr[(term['name'], db, term['goid'])] = distr.get((term['name'], db, term['goid']), 0) + len(dblink_ids)
                    else:
                        distr[(term['name'], db, term['goid'])] = distr.get((term['name'], db, term['goid']), 0) + 1

    remove_small_terms2(distr, small_num)

    return distr

##five statistical test functions using R
def hyper(q, m, n, k, is_tail = False):
    if m < q: return 0
    if is_tail:
        return robjects.r.phyper(q, m, n, k)[0]
    else:
        return 1 - robjects.r.phyper(q - 1, m, n, k)[0]

def binom_test(m1, n1, m2, n2, **kargs):
    p = float(m2) / n2
    alternative = kargs.get('alternative', 'greater')
    cmd = 'binom.test(c(%d, %d), p = %f, alternative = "%s")' \
          % (m1, n1 - m1, p, alternative)
    return robjects.r(cmd)[2][0]

def chisq_test(m1, n1, m2, n2, **kargs):
    ##m1         m2
    ##n1 - m2    n2 - m2
    cmd = 'chisq.test(matrix(c(%d, %d, %d, %d), nc = 2))' \
          % (m1, n1 - m1, m2, n2 - m2)
    res = robjects.r(cmd)
    if res[7][0] > 0:
        return res[2][0] / 2
    else:
        return 1 - res[2][0] / 2

def improved_chisq_test(m1, n1, m2, n2, **kargs):
    max_e = 10
    s1r, s2r, s1c, s2c, t = m1 + m2, (n1 + n2 - m1 - m2), n1, n2, float(n1 + n2)
    es = (s1r * s1c / t, s1r * s2c / t, s2r * s1c / t, s2r * s2c / t)
    for e in es:
        if e < max_e:    ##invoke Fisher's Exact Test when sample frequencies below 10
            return fisher_test(m1, n1, m2, n2, **kargs)
    return chisq_test(m1, n1, m2, n2, **kargs)

def fisher_test(m1, n1, m2, n2, **kargs):
    alternative = kargs.get('alternative', 'greater')
    cmd = 'fisher.test(matrix(c(%d, %d, %d, %d), nc = 2), alternative = "%s")' \
          % (m1, n1 - m1, m2, n2 - m2, alternative)
    return robjects.r(cmd)[0][0]

##three FDR correction functions
##Benjamini and Hochberg (1995)
def fdr_BH(ps):
    corrected_ps = list(robjects.r['p.adjust'](robjects.FloatVector(ps), method = 'BH'))
    return corrected_ps

##Benjamini and Yekutieli (2001)
def fdr_BY(ps):
    corrected_ps = list(robjects.r['p.adjust'](robjects.FloatVector(ps), method = 'BY'))
    return corrected_ps

##QVALUE
def fdr_QVALUE(ps):
    try:
        null = open('/dev/null', 'w')
        oldout = sys.stdout
        sys.stdout = null
        olderr = sys.stderr
        sys.stderr = null
        robjects.r.library('qvalue')
        corrected_ps = list(robjects.r['qvalue'](robjects.FloatVector(ps))[2])
    except:
        corrected_ps = ['None' for i in ps]
    finally:
        sys.stdout = oldout
        sys.stderr = olderr
    return corrected_ps

##two basic classes for test
class OneSampleTest:

    def __init__(self, *samples):
        self.samples = samples

    def calc_pvalue(self, *args):
        raise NotImplemented

    def __call__(self):
        raise NotImplemented

class TwoSampleTest(OneSampleTest):

    def __init__(self, *samples):
        if len(samples) != 2:
            raise TypeError, '%s expects 2 arguments, but gets %d' \
                % (self.__class__, len(*samples))
        OneSampleTest.__init__(self, *samples)

    def calc_pvalue(self, m1, n1, m2, n2):
        raise NotImplemented

    def __call__(self, result = None):
        sample1, sample2 = self.samples[:2]
        if isinstance(sample2.dis, Distr):
            for term in sample1.dis.keys():
                if not sample2.dis.has_key(term):
                    sample1.dis.pop(term)
            n1, n2 = sample1.gnum, sample2.gnum
            for term, queries in sample1.dis.items():
                m1, m2 = len(queries), len(sample2.dis[term])
                result.append(list(term) + ['%d' % m1, '%d' % m2, self.calc_pvalue(m1, n1, m2, n2), '|'.join(queries)])
        else:
            for term in sample1.dis.keys():
                if not term in sample2.dis.keys():
                    sample1.dis.pop(term)
            n1, n2 = sample1.gnum, sample2.gnum
             
            for term, queries in sample1.dis.items():
                m1, m2 = len(queries), sample2.dis[term]
                try:
                    result.append(list(term) + ['%d' % m1, '%d' % m2, self.calc_pvalue(m1, n1, m2, n2), '|'.join(queries)])
                except:
                    pass
        return result

##four test classes
class BinomTest(TwoSampleTest):

    def calc_pvalue(self, m1, n1, m2, n2, **kargs):
        return binom_test(m1, n1, m2, n2, **kargs)

class HyperTest(TwoSampleTest):

    def calc_pvalue(self, m1, n1, m2, n2, **kargs):
        return hyper(m1, m2, n2 - m2, n1, **kargs)

class ChiTest(TwoSampleTest):

    def calc_pvalue(self, m1, n1, m2, n2, **kargs):
        if m1 > m2:
            raise exception.StatError, 'Input is not a subset of background. Try binomial test.'
        return improved_chisq_test(m1, n1, m2, n2, **kargs)

class FisherTest(TwoSampleTest):

    def calc_pvalue(self, m1, n1, m2, n2, **kargs):
        if m1 > m2:
            raise exception.StatError, 'Input is not a subset of background. Try binomial test.'
        return fisher_test(m1, n1, m2, n2, **kargs)

class TestResult:

    def __init__(self, title = None):
        if title:
            self.title = title
        else:
            self.title = []
        self.result = []

    def __iter__(self):
        return iter(self.result)

    def __getitem__(self, i):
        return self.result[i]

    def __repr__(self):
        res = ''
        if self.title:
            res += '#%s\n' % '\t'.join(self.title)
        for feature in self.result:
            res += '\t'.join(map(str, feature)) + '\n'
        return res

    def __str__(self):
        return self.__repr__()

    def add_column(self, fieldname, fieldpos, column):
        assert(len(self.result) == len(column))
        self.title.insert(fieldpos, fieldname)
        for i in range(len(self.result)):
            self.result[i].insert(fieldpos, column[i])

    def add_fdr(self, method):
        do_fdr_func = globals()['fdr_%s' % method]
        corrected_ps = do_fdr_func([i[-2] for i in self.result])
        self.add_column('Corrected P-Value', -1, corrected_ps)

    def add_row(self, row):
        assert(len(row) == len(self.title))
        self.result.append(row)

    def add_title(self, title):
        self.title = title

    def append(self, row):
        self.add_row(row)

    def reverse(self):
        self.result.reverse()

    def sort(self, key = -1, order = 0):
        self.result.sort(list_cmp(key))
        if order:
            self.result.reverse()

class list_cmp:

    def __init__(self, i):
        self.i = i

    def cmp(self, x, y):
        '''select the criteria for comparison
        '''
        if x[self.i] < y[self.i]:
            return -1
        elif x[self.i] == y[self.i]:
            return 0
        else:
            return 1

    def __call__(self, x, y):
        return self.cmp(x, y)

