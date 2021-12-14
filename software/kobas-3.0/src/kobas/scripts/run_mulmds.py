# modified by Chen Ai, add log function for rna-seq data type
import os
from optparse import OptionParser
from kobas import config
kobasrc = config.getrc()

from kobas import snbea
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()


def config_option():
    usage = 'Usage: %prog -e expression_file -t datatype -p phenotype_file -o output_directory -s species [-d dbtype] [-i idtype] [-r normalize] [-m methods] [-n permuation_num] [-u min_size] [-l max_size] [-k kobas_home]'
    p = OptionParser(usage)
    p.add_option('-e', '--exp_file', dest = 'exp', action = 'store',type = 'string', 
    	help = 'Expression matrix. A tab delimited text file format that contains expression values. Columns correspond to samples,rows correspond to genes; Header=\'Genes\'\\t \\t sample names ,  Row names = gene ID. The type of the id should be assigned using -i .')
    p.add_option('-t', '--data_type', dest = 'datatype', action = 'store', type = 'string', 
    	help = 'Expression data type.  Use \'ma\' for microarray and \'rseq\' for RNA-seq data.')
    p.add_option('-r', '--normalize', dest = 'norm', action = 'store', default='n' ,type = 'string', 
        help = 'If you need normalization for expression data, \'-r y\', else \'-r n\' . Default \'n\', kobas supposes your data has not been normalized. Normalization is executed before differential expression, so this function is adapt to microarray data. Please note that normalization will not be done if your data type is RNA-Seq .')
    p.add_option('-p', '--phenotype_file',dest = 'phe',action = 'store',type = 'string', 
    	help = 'Phenotype data. One column file containing annotation information for the samples in expression matrix. The number of rows/samples in this file should match the number of columns/samples of the expression matrix. The column is reserved for a *BINARY* group assignment.  Use \'0\' and \'1\'for unaffected (controls) and affected (cases) sample class, respectively.')
    p.add_option('-o','--output_directory', dest='out',action = 'store',type = 'string',help = 'Directory of output files, must be existed, e.g. -o /path/to/output/')
    p.add_option('-s', '--species', dest = 'species', action='store', type = 'string',help = 'Input species abbreviation. For example: hsa for Homo sapiens, mmu for Mus musculus, dme for Drosophila melanogaster, ath for Arabidopsis thaliana, sce for Saccharomyces cerevisiae and eco for Escherichia coli K-12 MG1655, e.g. -s hsa')
    p.add_option('-d','--dbtype',dest = 'dbtype',action = 'store', default = 'K', type = 'string',
        help = 'Databases for selection, [K/R/B/p/k/G/S] 1-letter abbreviation, each time only run 1 kind of database. K for KEGG PATHWAY, R for Reactome, B for BioCyc, p for PANTHER, k for KEGG DISEASE, and G for Gene Ontology, S for GO slim. Default: K .')
    p.add_option('-i','--idtype',dest = 'idtype',action = 'store',default = 'ncbigene',
        type = 'string',help = 'Type of gene id in first column of expression matrix. Options: %s . Default: ncbigene. e.g. -i ncbigene' % (', '.join(DBLINKS.keys())))
    p.add_option('-m', '--methods', dest = 'methods', action = 'store' , default = 'GB/GS/GSA/P/GP/S/C/GA/PL', type='string',
            help='Methods to do enrichment analysis. Input 2-3 letter abbreviation separated by "/": GA for gage, GB for globaltest, GS for GSEA, GSA for gsa, P for padog, PL for plage, GP for GANPA, C for cepa, S for SAFE. Gene set-based methods: GA/GB/GS/GSA/P/PL/S; net-based methods: GP/C. Default: [GB/GS/GSA/P/GP/GE/C/GA/PL] .')
    p.add_option('-n','--npermutation',dest = 'nperm',action = 'store',default = 1000.0,
        type = 'float',help = 'Number of permutations of the expression matrix to estimate the null distribution. Defaults to 1000.')
    p.add_option('-u','--min_size',dest = 'min_size',action = 'store',default = 15,
                 type = 'int',help = 'Minimum size (in genes) for gene sets to be considered (default: 15)' )
    p.add_option('-l','--max_size',dest = 'max_size',action = 'store',default = 500,
                 type = 'int',help = 'Maximum size (in genes) for gene sets to be considered (default: 500)')
    p.add_option(
        '-k', '--kobashome', dest = 'kobas_home', default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to kobas_home, which is parent directory of sqlite3/ and seq_pep/ , default value is read from ~/.kobasrcwhere you set before running kobas. If you set this parameter, it means you set "kobasdb" and "blastdb" in this following directory. e.g. "-k /home/user/kobas/", means that you set kobasdb = /home/user/kobas/sqlite3/ '
    )
    opt,args=p.parse_args()
    return (p,opt,args)

DBLINKS = {'ncbigene': 'entrez_gene_id', 'ensembl': 'ensembl_gene_id'}
METHODS = {'GB':'globaltest',  'GS':'gsea', 'GSA':'gsa', 'P':'padog', 'GP':'ganpa', 'GE':'ggea', 'C':'cepa', 'GA':'gage', 'PL':'plage', 'O':'ora','S':'safe'}

#SBEA = ("ora","globaltest","samgs","gsea","padog","gsa")
#NBEA = ("ganpa","ggea","cepa","pathnet", "spia", "topologygsa")
#OTHERS = ('gage','plage')

svc_ma_mds = ('cepa', 'gage', 'ganpa', 'globaltest', 'gsa', 'gsea','padog', 'plage', 'safe')
svc_rseq_mds =  ('cepa', 'gage', 'ganpa', 'ggea', 'globaltest','plage')

opt_parser, opt, args = config_option()

if opt.kobas_home:
    kobasrc['kobas_home'] = opt.kobas_home
    kobasrc['kobasdb'] = opt.kobas_home + '/sqlite3/'
    kobasrc['gmt'] = opt.kobas_home + '/gmt/'
    kobasrc['grn'] = opt.kobas_home + '/grn/'
    kobasrc['model'] = opt.kobas_home + '/model/'


if not opt.out:
    opt_parser.error('Option -o must be assigned.\n')


if opt.species and opt.dbtype and opt.idtype:
    if len(opt.dbtype)>1:
        opt_parser.error('Option -d database only allow 1 kind database.')
    gmtfiles = [kobasrc['kobas_home'] +'/gmt/' + ('.').join([opt.species, opt.dbtype, DBLINKS[opt.idtype],'gmt'])]
#    dbs = opt.dbtype.split('/')
#    gmtfiles = []
#    for db in dbs:
#    	gmtfile = kobasrc['kobas_home'] +'/gmt/' + ('.').join([opt.species, db, DBLINKS[opt.idtype],'gmt'])
#    	gmtfiles.append(gmtfile)
else :
    opt_parser.error('Option -s -d -i must be assigned. \n')

if not opt.datatype in ('ma','rseq'):
    opt_parser.error('Only support data come from microarray (\'ma\') and RNA-Seq (\'rseq\'). Please check -t ')

r = robjects.r
biobase = importr("Biobase")
enbrowser = importr("EnrichmentBrowser")

if opt.methods:
    tmp_md = opt.methods.split('/')
    methods = [ METHODS[x] for x in tmp_md ]
    #if opt.datatype == 'rseq':
    #    methods = tuple(set(methods) & set(svc_rseq_mds))
    #    print "Warning: for RNA-Seq datatype, only these methods are supported: %s " % ' '.join(svc_rseq_mds)
    sbea_mds, nbea_mds, other_mds = snbea.select_methods(methods)
if nbea_mds :
    grnfile = kobasrc['kobas_home'] +'/grn/' + ('.').join([opt.species, DBLINKS[opt.idtype], 'grn'])
    grn = snbea.get_grn(grnfile)
#    grn = pd.read_table(grnfile, header=None, names=["FROM","TO","TYPE"], dtype=str)
#    grn = pandas2ri.py2ri(grn)
#    grn = r['data.frame'](grn, stringsAsFactors=False)
#    grn = r['as.matrix'](grn)

gmt_dict = {} # each values is a pandas dataframe contain description of gene sets
gs_dict = {}  # each value is a r list contain genes of gene sets
for gmt in gmtfiles:
    gset_des, gset_genes = snbea.read_gmt_file(gmt)
    gmt_dict[gmt] = pd.DataFrame({'GENE.SET':gset_des.keys(), 'TITLE':gset_des.values()})
    genes = [ robjects.StrVector(tuple(x)) for x in gset_genes.values() ]
    gs_dict[gmt] = pandas2ri.ListVector(dict(zip(gset_genes.keys(), genes)))     
    gmt_dict[gmt], gs_dict[gmt] = snbea.filter_gsets(gmt_dict[gmt], gs_dict[gmt], opt.min_size, opt.max_size) 

gset_df = pd.concat(gmt_dict.values())
if len(gs_dict) >1 :
    gsets = r['append'](*gs_dict.values())
else:
    gsets = gs_dict.values()[0]

expdata = snbea.read_exp_file(exp=opt.exp)
pdata = snbea.get_pdata(phe=opt.phe, expdata=expdata)
print 'Normalize and Do differently expressed gene analysis... \n'
exset = snbea.get_exset(expdata,pdata, opt.datatype, opt.norm)
print 'Output the differently expressed genes file named as \'degenes.tsv\' \n'
snbea.output_de_genes(exset, opt.out)


res = {}
tbl = {}

for sbea_md in sbea_mds:
    try:
        print '#'*20
        print ' '*5+sbea_md
        print '#'*20
        res[sbea_md] = enbrowser.sbea(method=sbea_md, eset=exset, gs=gsets, alpha=0.999, perm=opt.nperm, out_dir=opt.out)
        tbl[sbea_md] = snbea.enres_add_annotation(res[sbea_md], gset_df)
        tbl[sbea_md] = tbl[sbea_md].sort_values(by=['P.VALUE'])
        filenm = opt.out + sbea_md+ '.tsv'
        tbl[sbea_md].to_csv(filenm, sep="\t", header=True, index=False)
    except:
        print 'Error: '+ sbea_md +' have not successfully run.'
        pass

for nbea_md in nbea_mds:
    try :        
        print '#'*20
        print ' '*5+nbea_md
        print '#'*20
        res[nbea_md] = enbrowser.nbea(method=nbea_md, eset=exset, gs=gsets, alpha=0.9999, grn=grn, out_dir=opt.out, perm=opt.nperm)
        tbl[nbea_md] = snbea.enres_add_annotation(res[nbea_md], gset_df)
        tbl[nbea_md] = tbl[nbea_md].sort_values(by=['P.VALUE'])
        filenm = opt.out + nbea_md+ '.tsv'
        tbl[nbea_md].to_csv(filenm, sep="\t", header=True, index=False)
    except:
        print 'Error: '+ nbea_md + ' have not successfully run.'
        pass        


for other_md in other_mds:
    try:
        print '#'*20
        print ' '*5+ other_md
        print '#'*20
        if other_md == 'gage':
            tbl[other_md] = snbea.run_gage(exset, gsets, gset_df)
        if other_md == 'plage':
            tbl[other_md] = snbea.run_plage(exset, gsets, gset_df)
        filenm = opt.out + other_md + '.tsv'
        tbl[other_md].to_csv(filenm, sep="\t", header=True, index=False)
    except:
        print 'Error: '+ other_md + ' have not successfully run.'
        pass


######## combine all the result with rank of SVM classifier ########

from sklearn.externals import joblib
#if opt.datatype == 'ma':
#    svc_mds = svc_ma_mds
#    svc = joblib.load(kobasrc['model'] + '/svc.for_ma.pkl')
#else :
#    svc_mds = svc_rseq_mds
#    svc = joblib.load(kobasrc['model'] + '/svc.for_rseq.pkl')

svc_features = [
    'cepa',
    'gage',
    'ganpa',
    'globaltest',
    'gsa',
    'gsea',
    'padog',
    'plage',
    'safe']   #### this order can't be change, each SVC dimension should be corresponding

svc_mds = svc_ma_mds
svc = joblib.load(kobasrc['model'] + '/svc.pkl')
md_rank={}
rank_tb = gset_df
if set(tbl.keys()) >= set(svc_mds):
    for md in svc_mds:
        if md == 'gage':
            pval = 'p.val'
        elif md == 'plage':
            pval = 'P.Value'
        else:
            pval = 'P.VALUE'
        tbl[md]['rank'] = tbl[md][pval].rank(method = 'min', ascending=True)
        md_rank[md] = pd.DataFrame({'GENE_SET': tbl[md]['GENE.SET'], md:tbl[md]['rank']})        
        rank_tb = rank_tb.merge(md_rank[md], left_on = 'GENE.SET', right_on = 'GENE_SET', how = 'outer')       
        rank_tb = rank_tb.drop('GENE_SET',axis=1)
    rank_arr = rank_tb.iloc[:,2:]
    rank_arr = rank_arr[svc_features]   ### column order follow as svm's features 
    rank_arr = rank_arr /  (1.0*rank_arr.shape[0]) ### divide the total of the gene sets invert to rank percent
    rank_arr = rank_arr.fillna(1.0)
     
    # start classify
    comb_cls = svc.predict(rank_arr).astype(int)
    comb_prb = svc.predict_proba(rank_arr)   # 2 class, so 2-d array : probability to predict as 0 and 1 
    comb_prb = np.max(comb_prb,axis=1)   #  the probability to predict correctly
    comb_dis = svc.decision_function(rank_arr)
    rank_out = rank_tb.iloc[:,:2]    
    tmp = pd.DataFrame(np.column_stack([comb_cls,comb_prb,comb_dis]))
    rank_out['ENT'] = tmp.iloc[:,0]==1
    rank_out = pd.concat([rank_out, tmp.iloc[:,1:]], axis=1, ignore_index=True)
    rank_out.columns = ['GENE_SET','NAME','ENRICHMENT_RES','PROBABILITY','ENRICH_SCORE']
    rank_out = rank_out.sort_values('ENRICH_SCORE', ascending=False)
    rank_out.to_csv(opt.out +'combination_results.tsv',  sep="\t", header=True, index=False)

else:
    print 'Warning: Sorry, it is not able to do multiple methods results combination with these input methods. All of these methods results are required. ' + ', '.join(svc_mds)

