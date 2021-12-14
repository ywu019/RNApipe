from kobas import config
import numpy as np
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

from kobas import exception


def select_methods(mds):
    sbea_mds = []
    nbea_mds = []
    other_mds = []
    for md in mds:
    	if md in all_sbea_mds:
	    sbea_mds.append(md)
	elif md in all_nbea_mds:
            nbea_mds.append(md)
        else:
            other_mds.append(md)
    return sbea_mds, nbea_mds, other_mds


r = robjects.r
enbrowser = importr("EnrichmentBrowser")
biobase = importr("Biobase")
edgeR = importr("edgeR")

all_sbea_mds = tuple(enbrowser.sbea_methods())
all_nbea_mds = tuple(enbrowser.nbea_methods())

def read_gmt_file(gmt):
    gset_des={}    #description of gene sets
    gset_genes = {}  # genes of gene sets
    gmt_file = open(gmt,'r')
    line = gmt_file.readline()
    try :
        while line :
            if '#' in line:
            	line = gmt_file.readline()
            	continue
            words = line.rstrip('\n').split('\t')

            if len(words)<3 :
                raise exception.GmtformatErr, "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file, a tab delimited file are allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats "

            if words[0] in gset_des.keys():
                continue
                print "Warning Message from KOBAS: \n These are at leaset two gene sets share the same gene set ID. Please check. We only use the gene set first emerges."
        
            gset_des[words[0]] = words[1]
            gset_genes[words[0]] = words[2:]
            line = gmt_file.readline()
        gmt_file.close()
    except IOError,e:
        print "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file are allowed."
        gmt_file.close()
        sys.exit(e)
    return gset_des, gset_genes


def read_exp_file(exp):
    # the first column must be gene name!
    r('''
    read_exp <- function(exp){
      expdf = read.delim(exp, stringsAsFactors=F, header=T)
      expdata = as.matrix(expdf[,2:ncol(expdf)])
      rownames(expdata) <- expdf[,1]
      return(expdata)
    }
    '''
    )
    read_exp = r['read_exp']
    return read_exp(exp)


def get_pdata(phe,expdata):
    with open(phe,'r') as f:
        lines = f.readlines()
        pheno = [x.rstrip() for x in lines]
    pheno = robjects.IntVector(pheno)
    samples = r['colnames'](expdata)
    pdata = r['data.frame'](SAMPLE= samples, GROUP=pheno, stringsAsFactors=False)
    pdata.rownames = samples
    pdata = biobase.AnnotatedDataFrame(pdata)
    return pdata


def filter_gsets(gmt, gs, min_size, max_size):
    gset_size = r['lapply'](gs ,r.length)
    ltlist = list(gset_size.ro > min_size)
    mtlist = list(gset_size.ro < max_size)
    func = lambda x,y : True if x and y else False
    adjpy = map(func, ltlist, mtlist)
    adjr  = robjects.BoolVector(adjpy)
    gmt = gmt[adjpy]
    gs = gs.rx(adjr)
    return gmt,gs

def get_grn(grnfile):
    grn_py = pd.read_table(grnfile, header=None, names=["FROM","TO","TYPE"], dtype=str)
    grn = pandas2ri.py2ri(grn_py)
    grn = r['data.frame'](grn, stringsAsFactors=False)
    grn = r['as.matrix'](grn)
    return grn

def get_exset(expdata,pdata,datype,norm):
    ### differently express genes #######
    if datype=='ma':
        print "-"*20
        print "Microarray data"
        print "-"*20
        exset = biobase.ExpressionSet(assayData=expdata, phenoData=pdata)
        # normalize:
        if norm =='y':
            print "Doing Normalization"
            exset = enbrowser.normalize(exset, data_type=datype)
    if datype=='rseq':
        print "-"*20
        print "RNA-SEQ"
        print "-"*20
        if norm != 'y':
            expdata = r['log2'](expdata.ro + 1.0)
        else:
            expdata = edgeR.cpm(expdata, log=True, prior_count=3)
        exset = biobase.ExpressionSet(assayData=expdata, phenoData=pdata)    
    exset = enbrowser.de_ana(exset, de_method = 'limma')   # it will check datatype auto, decimal: microarray, integer:RNA-Seq;(So, after log tranformation, RNA-Seq is regarded as ma). For RNA-Seq read count, TMM and VOOM is done, see EnrichmentBrowser source code deAna.R
    return exset

def run_plage(exset, gsets, gset_df):
    plagepkg = importr("GSVA")
    limma = importr("limma")
    phdata = biobase.pData(exset)
    phdata_py = pandas2ri.ri2py_dataframe(phdata)
    group_df = pd.DataFrame({'ctrl':[1 for i in range(len(phdata_py)) ], 'case':phdata_py.iloc[:,1]})
    group_df = pandas2ri.py2ri_pandasdataframe(group_df)
    design = r['as.matrix'](group_df)
    plage_res = plagepkg.gsva(exset, gsets, method="plage", rnaseq=False)[0]
    fit = limma.lmFit(plage_res, design)
    fit = limma.eBayes(fit)
    plage_tbl = limma.toptable(fit, coef=2, number=len(gsets))
    plage_tbl = pandas2ri.ri2py_dataframe(plage_tbl)
    plage_tbl['GENE.SET'] = plage_tbl.index
    plage_tbl = pd.merge(gset_df, plage_tbl, how='inner', sort=False)
    plage_tbl = plage_tbl.sort_values(by=['P.Value'])
    return plage_tbl

def run_gage(exset, gsets, gset_df):
    gagepkg = importr("gage")
    expdata = biobase.exprs(exset)
    phdata = biobase.pData(exset)
    phdata = phdata.rx(2)
    phdata_py = pandas2ri.ri2py_dataframe(phdata)
    phdata_py.index = range(len(phdata_py))
    phdata_py.index = phdata_py.index + 1
    ctrl = phdata_py.iloc[:,0] == 0
    ctrl_num = pandas2ri.IntVector(ctrl[ctrl].index)
    case_num = pandas2ri.IntVector(ctrl[~ctrl].index)
    if len(ctrl_num) == len(case_num):
        gage_res = gagepkg.gage(exprs=expdata, gsets = gsets, ref=ctrl_num, samp=case_num, same_dir=False)
    else:
        gage_res = gagepkg.gage(exprs=expdata, gsets = gsets, ref=ctrl_num, samp=case_num, same_dir=False, compare='unpaired')
    gage_res_dict = dict(zip(gage_res.names, list(gage_res))) 
    gage_res_rdf = gage_res_dict['greater']
    gage_res_rdf = r['data.frame'](gage_res_rdf)
    gage_res_df = pandas2ri.ri2py_dataframe(gage_res_rdf)
    gage_res_df['GENE.SET']=gage_res_df.index
    gage_final = pd.merge(gset_df,gage_res_df,how='inner',sort=False)
    gage_final = gage_final.sort_values(by=['p.val'])
    return gage_final


##### output de genes	
def output_de_genes(exset, outpath):
    fdata = biobase.fData(exset)
    fdata = pandas2ri.ri2py_dataframe(fdata)
    fdata = fdata.sort_values(by=['ADJ.PVAL'])
    fdata.index.name = 'GENE'
    fdata.to_csv(outpath +'degenes.tsv',sep='\t')


def enres_add_annotation(res, gset_df) :
    res_dict = dict(zip(res.names, list(res)))
    tbl = res_dict['res.tbl']
    tbl = r['data.frame'](tbl, stringsAsFactors=False)
### add title to pathways
    tbl_df = pandas2ri.ri2py_dataframe(tbl)
    tbl_df = pd.merge(gset_df, tbl_df, how='inner',sort=False)
    return tbl_df


