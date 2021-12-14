#!/usr/bin/env python
import sys,time
from optparse import OptionParser
from kobas import gsea,config,dbutils,exception

def config_option():
    usage = 'Usage: %prog -e expression_file -c cls_file -o output_directory [-g gmt_file] [-s species] [-d dbtype] [-i idtype] [-n npermutation] [-w weighted_score_type] [-m method] [-u min_size] [-l max_size] [-k kobas_home] [-q kobasdb]'
    p = OptionParser(usage)
    p.add_option('-e','--exp_file',dest = 'exp',action = 'store',type = 'string',help = 'input gct file or txt file, e.g. -e abc.gct')
    p.add_option('-c','--cls_file',dest = 'cls',action = 'store',type = 'string',help = 'input cls file, e.g. -c abc.cls')
    p.add_option('-o','--output_directory',dest = 'out',action = 'store',type = 'string',help = 'directory of output files, must be existed,e.g. -o /path/to/output')
    p.add_option('-g','--gmt_file',dest = 'gmt',action = 'store',type = 'string',help = 'input gmt file if without species, e.g. -g 123.gmt')
    p.add_option('-s','--species',dest = 'species',action = 'store',type = 'string',help = 'input species abbreviation if without gmt file,  (for example: ko for KEGG Orthology, hsa for Homo sapiens, mmu for Mus musculus, dme for Drosophila melanogaster, ath for Arabidopsis thaliana, sce for Saccharomyces cerevisiae and eco for Escherichia coli K-12 MG1655), e.g. -s hsa')
    p.add_option('-d','--dbtype',dest = 'dbtype',action = 'store', default = 'P:K/n/b/R/B/p', type = 'string',
                 help = 'databases for selection, Database category: Database type; Database category(P, for pathway; D, for disease when -s hsa; G for Gene Ontology), Database type: 1-letter abbreviation separated by "/": K for KEGG PATHWAY, n for PID, b for BioCarta, R for Reactome, B for BioCyc, p for PANTHER, o for OMIM, k for KEGG DISEASE, f for FunDO, and G for Gene Ontology. Only one database category is allowed. [P:K/n/b/R/B/p] [D:o/k/f] [G] Default: P:K/n/b/R/B/p, e.g. D:k/o')
    p.add_option('-i','--idtype',dest = 'idtype',action = 'store',default = 'id:ncbigene',
                 type = 'string',help = 'choose idtype if species but not gmt file is assigned. Options: %s , default: id:ncbigene. e.g. -i id:uniprot' % (', '.join(DBLINKS.keys())))
    p.add_option('-n','--npermutation',dest = 'nperm',action = 'store',default = 1000,
                 type = 'int',help = 'Times of random permutations to perform')
    p.add_option('-w','--weighted_score_type',dest = 'wst',action = 'store',default = 1.0,
                 type = 'float',help = 'Type of score: weight. Options: 0 (unweighted = Kolmogorov-Smirnov), 1.0 (weighted), and 2.0 (over-weighted),or other numbers between 0-2, e.g. -w 1.5 (default 1.0)')
    p.add_option('-m','--method',dest= 'method',action = 'store',default = 'snr',
                 type = 'string',  help = 'Method type used in ranking genes in list. Options: snr (signal to noise calculation, uses the difference of means scaled by the standard deviation), ttest (tTest, uses the difference of means scaled by the standard deviation and number of samples), default : snr , e.g. -m ttest')
    p.add_option('-u','--min_size',dest = 'min_size',action = 'store',default = 15,
                 type = 'int',help = 'Minimum size (in genes) for database gene sets to be considered (default: 15)' )
    p.add_option('-l','--max_size',dest = 'max_size',action = 'store',default = 500,
                 type = 'int',help = 'Maximum size (in genes) for database gene sets to be considered (default: 500)')
## add to replace kobasrc
    p.add_option(
        '-k', '--kobashome', dest = 'kobas_home', default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to kobas_home, which is parent directory of sqlite3/ and seq_pep/ , default value is read from ~/.kobasrcwhere you set before running kobas. If you set this parameter, it means you set "kobasdb" and "blastdb" in this following directory. e.g. "-k /home/user/kobas/", means that you set kobasdb = /home/user/kobas/sqlite3/ '
    )
    p.add_option(
        '-q','--kobasdb',dest = 'kobasdb',default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to sqlite3/, default value is read from ~/.kobasrc where you set before running kobas, e.g. "-q /kobas_home/sqlite3/"'
    )    
    opt,args=p.parse_args()
    return (p,opt,args)

DBLINKS = {'id:ncbigene': 'entrez_gene_id', 'id:ncbigi': 'gi', 'id:uniprot': 'uniprotkb_ac', 'id:ensembl': 'ensembl_gene_id'}
METHODS = ('snr','ttest')

def verify_exp_file(exp_file, sample_num):
    if '.txt' in exp_file :
        txt = exp_file
        (gene_list,gene_num,expr_data,title) = gsea.read_txt_file(txt=txt, sample_num = sample_num)
    elif '.gct' in exp_file :
        gct = exp_file
        (gene_list,gene_num,expr_data,title) = gsea.read_gct_file(gct=gct, sample_num = sample_num)
    else:
        raise exception.FormatErr, "Error Message: \nWrong format. Only TXT and GCT format is allowed, which is same as GSEA."
    return gene_list,gene_num,expr_data,title 

opt_parser, opt, args = config_option()
kobasrc = config.getrc()
if opt.kobas_home:
    kobasrc['kobas_home'] = opt.kobas_home
    kobasrc['kobasdb'] = opt.kobas_home + '/sqlite3/'
if opt.kobasdb:
    kobasrc['kobasdb'] = opt.kobasdb +'/'

if not opt.out:
    opt_parser.error('Option -o must be assigned.\n')    
### get constants and verify
if opt.species:
    speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + opt.species + '.db')


if not opt.method in METHODS:
    opt_parser.error('This method is not allowed when running GSEA, please use -h to read options help message')


if opt.idtype :
    idtype = DBLINKS[opt.idtype]
else:
    idtype = DBLINKS['id:ncbigene']

### start to read files/databases
if opt.cls:
    (labels,phnameA,phnameB,sample_num,sample0,sample1) =gsea.read_cls_file(cls=opt.cls)
else :
    opt_parser.error('Option -c must be assigned.\n')

if opt.exp :
    (gene_list,gene_num,expr_data,title) = verify_exp_file(exp_file=opt.exp, sample_num = sample_num)   
else :
    opt_parser.error('Option -e must be assigned.\n')

if opt.gmt:
    (gset_name,gset_des,gset_genes) = gsea.read_gmt_file(gmt=opt.gmt)
    (hit_matrix_filtered, hit_genes_filterd, hit_sum_filtered, gset_name_filtered, gset_des_filtered) = \
        gsea.get_hit_matrix(gene_list, gene_num, gset_name, gset_des,gset_genes, opt.min_size, opt.max_size)
    all_gset_num = len(gset_name)
elif opt.species:
    pre_gene_list, pre_gene_num, pre_expr_data = gene_list, gene_num, expr_data
    (gene_list, gene_num, expr_data) = gsea.pick_genes(pre_gene_list,pre_gene_num,pre_expr_data, speciesdb, idtype, opt.out, title)
    (all_gset_num, hit_matrix_filtered, hit_genes_filterd, hit_sum_filtered, gset_name_filtered, gset_des_filtered) = \
            gsea.read_gmt_db(gene_list,speciesdb,opt.dbtype,idtype,opt.min_size,opt.max_size, opt.species)
else:
    opt_parser.error('Either option -g or option -s must be assigned')


(sort_r, sort_gene_index) = gsea.rank_pro(lb=labels,md=opt.method,expr_data=expr_data,sample0=sample0,sample1=sample1)
(ES,ESidx,RES) = gsea.ES_all(sort_r=sort_r, sort_gene_index=sort_gene_index,hit_matrix_filtered=hit_matrix_filtered, weighted_score_type=opt.wst, gene_num = gene_num)
es_null = gsea.ES_null(lb=labels, times=opt.nperm, method=opt.method, sample0=sample0, sample1=sample1, hit_matrix_filtered = hit_matrix_filtered,
                       weighted_score_type = opt.wst,expr_data=expr_data,gene_num=gene_num)
pval = gsea.nominal_p(es = ES, es_null = es_null)
(NES,nes_null) = gsea.normalized(ES,es_null)
fdr = gsea.fdr_cal(nes_obs=NES, nes_null=nes_null)
(type0,type1) = gsea.output_set(result_path=opt.out, phnameA=phnameA,phnameB=phnameB,gset_name_filtered=gset_name_filtered,
                                gset_des_filtered=gset_des_filtered,hit_matrix_filtered=hit_matrix_filtered,ES=ES, NES=NES,nompval=pval, FDR=fdr)

gsea.gene_info(result_path=opt.out, gene_list=gene_list, phnameA=phnameA, phnameB=phnameB, type0=type0, type1=type1,
               hit_matrix_filtered=hit_matrix_filtered, sort_r=sort_r, sort_gene_index=sort_gene_index, es_idx=ESidx, RES=RES )

gsea.plots_RES(result_path=opt.out , phnameA=phnameA, phnameB=phnameB, type0=type0,type1=type1,RES=RES, es_idx=ESidx,
               hit_matrix_filtered=hit_matrix_filtered,sort_gene_index=sort_gene_index, sort_r=sort_r)


method_dict = {'snr':'signal to noise ratio', 'ttest':'T-test'}

summary = open(opt.out + '/summary', 'w')
summary.write("#Summary of GSEA \n")
summary.write("#Input Gene number: " + str(pre_gene_num) +'; ' + "Analysis Gene Number: " + str(gene_num) + ' ; We only analyze genes that are matched with KEGG GENE. To see which genes are analysis, please read file \"Genes_Used_In_Analysis.txt\". To see which genes are not matched with KEGG GENE, please read the file \"Genes_abandoned.txt\". '+'\n' )
summary.write("#All Gene Sets :" + repr(all_gset_num) +"\n")
summary.write("#Gene set size filters (min="+ repr(opt.min_size) +", max=" +repr(opt.max_size) + ") resulted in filtering out " + repr(all_gset_num-len(ES))+ " gene sets. \n")
summary.write('#Times of random permutation on samples : '+ str(opt.nperm) +'\n')
summary.write('#The weight of score : ' + str(opt.wst) + '\n')
summary.write('#The method used in ranking genes in list : '+ method_dict[opt.method] +'\n\n')

summary.write("#The remaining " + repr(len(ES))+ " gene sets were used in the analysis.The number of significant sets are as follows:\n")
summary.write("Phenotype\tSamples\tEnrichment Set\tNominal P-value < 0.05\tFDR < 0.25\tGenes(Features)\n")
summary.write(phnameA + "\t" + repr(sample0) + "\t" + repr(len(type0)) + "\t" + repr(sum(type0['Nominal p-value']<0.05)) +"\t" + repr(sum(type0['FDR']<0.25)) +"\t"+repr(sum(sort_r>0))+"\n")
summary.write(phnameB + "\t" + repr(sample1) + "\t" + repr(len(type1)) + "\t" + repr(sum(type1['Nominal p-value']<0.05)) +"\t" + repr(sum(type1['FDR']<0.25)) +"\t"+repr(sum(sort_r<0))+"\n")
summary.write("Total\t" + repr(sample_num)+"\t"+repr(len(ES)) +"\t"+ repr(sum(type0['Nominal p-value']<0.05)+sum(type1['Nominal p-value']<0.05)) + "\t" + repr(sum(type0['FDR']<0.25)+sum(type1['FDR']<0.25)) + "\t" +repr(sum(sort_r>0)+sum(sort_r<0)) +"\n")
