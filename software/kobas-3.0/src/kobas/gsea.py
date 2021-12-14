#!/usr/bin/env python
import os,sys,math,warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot,savefig
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages


from kobas import config,dbutils,exception



def read_gmt_file(gmt):
    gset_name=[]    #name of gene sets
    gset_des=[]     #description of gene sets
    gset_genes=[]       #genes in gene sets
    gmt_file = open(gmt,'r')
    line = gmt_file.readline()
    try :
        while line :
            words = line.rstrip('\n').split('\t')
            if len(words)<3 :
                raise exception.GmtformatErr, "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file, a tab delimited file are allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats "

            if words[0] in gset_name:
                raise exception.GmtformatErr, "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file are allowed.\nDuplicate gene set names are not allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats "

            gset_name.append(words[0])
            gset_des.append(words[1])
            gset_genes.append(tuple(words[2:],))
            line = gmt_file.readline()
        gmt_file.close()
    except IOError,e:
        print "Error Message:\nBad gene set database format. Please check your gene set database format, only gmt file are allowed."
        gmt_file.close()
        sys.exit(e)
    gset_name = np.array(tuple(gset_name))
    gset_des = np.array(gset_des)
    gset_genes = np.array(gset_genes)

    return gset_name,gset_des,gset_genes


def get_hit_matrix(gene_list, gene_num, gset_name, gset_des,gset_genes, min_size, max_size):
    '''
        This function get a matrix of all gene sets
    '''

    try:
        hit_m =[]
        hit_genes=[]
        for i in range(len(gset_name)):
            hit_g = []
            is_hit = 0 # 0 -- not hit, 1 -- hit
            hit_array = []  # save  which position of sorted list is hitted
            for j in range(gene_num):
                gene = gene_list[j]
                #multi_gene = gene.split(' /// ')  # for some gene like 'KRT6E /// KRT6A /// KRT6C''
                #for each in multi_gene :
                if gene in gset_genes[i] :
                    is_hit = 1
                    hit_g.append(gene)
                hit_array.append(is_hit)
                is_hit = 0
            hit_genes.append(hit_g)
            hit_m.append(hit_array)
        hit_matrix = np.array(hit_m)
        if not 1 in hit_matrix:
            raise exception.GmtformatErr,"Error Message:\nNone of genes in input file(gct file) was matched with genes in gene sets. \nPlease check the gct file(-e), the gmt file(-g), the idtype(-i), the database type(-d)."

        hitsum = hit_matrix.sum(1)
        delindex = np.where((hitsum<min_size) | (hitsum>max_size))[0]
        hit_matrix_filtered = np.delete(hit_matrix,delindex,axis=0)
        hit_genes_filterd = np.delete(hit_genes,delindex,axis=0)
        hit_sum_filtered = np.delete(hitsum, delindex, axis=0)
        gset_name_filtered = np.delete(np.array(gset_name),delindex,axis=0)
        gset_des_filtered = np.delete(gset_des,delindex,axis=0)
        if len(hit_matrix_filtered)==0:
            raise exception.GmtvalueErr,"Error Message:\nAll gene sets"+repr(len(hit_matrix))+"have been filtered. \nPlease check the threshold and ceil of gene set size (values of min_size and max_size). "
    except ValueError,e:
        sys.exit(e)
    return hit_matrix_filtered, hit_genes_filterd, hit_sum_filtered, gset_name_filtered, gset_des_filtered

def read_gmt_db(gene_list,speciesdb, dbtype, idtype, min_size, max_size, species):

    gset_name=[]    #name of gene sets
    gset_des=[]     #description of gene sets
    hit_genes = []      #genes in gene sets none of use, only use genes both in list and gene set
    D = {'o': 'OMIM', 'k': 'KEGG DISEASE', 'f': 'FunDO'}  #don't support GAD and N, because lack of did

    try:
        dbcatg = dbtype.split(':')[0]
        if not dbcatg in ['P','D','G']:
            raise exception.GmtdatabaseErr,"Error Message: This database category you input is not allowed. Please choose 'P','D','G' for database category."

        records =[]
        gset_num = 0

        if dbcatg=='P':
            input_dbs = set(dbtype.split(':')[1].split('/'))
            avail_dbs = {}
            for database in speciesdb.pathwaydbs_from_abbr(species):
                avail_dbs[database[0]] = database[1]
            dbs = input_dbs.intersection(set(avail_dbs.keys()))
            if dbs:
                print  "Databases: %s"  % ', '.join([avail_dbs[db] for db in dbs])
            else:
                raise exception.GmtdatabaseErr, "No supported databases are selected. Supported databases are %s, but your input databases are: %s." % \
                ('/'.join(avail_dbs.keys()), dbtype)


            for db in dbs:
                tmp = speciesdb.allpathways(db).fetchall()
                num = speciesdb.pathwaynums(db)
                records += tmp
                gset_num += num
            if gset_num ==0 :
                raise exception.GmtdatabaseErr,"Error Message: \nFail to get information from gene set database."

            i = 0
            hit_matrix = np.array([[0.0 for s in range(len(gene_list))] for m in range(gset_num)])
            for record in records:
                hgene = []
                gset_name.append(record[1])     #id
                gset_des.append(record[2])      #name
                gset_gids = speciesdb.genes_from_pid(record[0])   #from table GenePathways
                for gid in gset_gids:
                    #print gid
                    #dblink_id = speciesdb.dblink_ids_from_gid(gid[0], idtype).fetchone()[0]
                    for each in speciesdb.dblink_ids_from_gid(gid[0], idtype):
                        dblink_id = each[0]
                    if dblink_id in gene_list:
                        hgene.append(dblink_id)
                        j = gene_list.index(dblink_id)
                        hit_matrix[i,j] = 1
                hit_genes.append(hgene)
                i+=1
            if not 1 in hit_matrix:
                raise exception.GmtdatabaseErr,"Error Message:\nNone of genes in input file(gct file) was matched with genes in gene sets. \nPlease check the gct file(-e), the species(-s), the idtype(-i), the database type(-d)."
        
        elif dbcatg=='D':
            input_dbs = set(dbtype.split(':')[1].split('/'))
            avail_dbs = {}
            if species != 'hsa':
                raise exception.GmtdatabaseErr,"Error Message: Disease is only supported for Homo Sapiens(-s hsa), not supported for this species. Please choose another database category, e.g. 'P','G'"
            dbs = input_dbs.intersection(set(D.keys()))

           # if dbs:
           #     print  'Databases: %s'  % ', '.join([avail_dbs[db] for db in dbs])
           # else:
            if not dbs:
                raise exception.GmtdatabaseErr, 'No supported databases are selected. Supported databases are %s, but your input databases are: %s.' % \
                ('/'.join(avail_dbs.keys()), dbtype)
            
            for db in dbs:
                tmp = speciesdb.alldiseases(db).fetchall()
                num = speciesdb.diseasenums(db)
                records += tmp
                gset_num += num
            if gset_num == 0:
                raise exception.GmtdatabaseErr,"Error Message: \nFail to get information from gene set database."

            i = 0
            hit_matrix = np.array([[0.0 for s in range(len(gene_list))] for m in range(gset_num)])
            for record in records:
                hgene = []
                #gset_name.append(record[1])     #id
                if record[1] :
                    gset_name.append(record[1])
                else:
                    #gset_name.append(record[2]) 
                    continue
                gset_des.append(record[2])      #name
                gset_gids = speciesdb.genes_from_did(record[0])   #from table GenePathways
                for gid in gset_gids:
                    #dblink_id = speciesdb.dblink_ids_from_gid(gid[0], idtype).fetchone()[0]
                    for each in speciesdb.dblink_ids_from_gid(gid[0], idtype):
                        dblink_id = each[0]
                    if dblink_id in gene_list:
                        hgene.append(dblink_id)
                        j = gene_list.index(dblink_id)
                        hit_matrix[i,j] = 1
                hit_genes.append(hgene)
                i+=1
            if not 1 in hit_matrix:
                raise exception.GmtdatabaseErr,"Error Message:\nNone of genes in input file(gct file) was matched with genes in gene sets. \nPlease check the gct file(-e), the species(-s), the idtype(-i), the database type(-d)."

        else :
            records = speciesdb.allgoterms().fetchall()
            gset_num = speciesdb.gotermnums()
            i=0
            hit_matrix = np.array([[0.0 for s in range(len(gene_list))] for m in range(gset_num)])
            for record in records:
                hgene = []
                gset_name.append(record[0])     #id
                gset_des.append(record[1])      #name
                gset_gids = speciesdb.genes_from_goid(record[0])   #from table GenePathways
                for gid in gset_gids:
                    #dblink_id = speciesdb.dblink_ids_from_gid(gid[0], idtype).fetchone()[0]
                    dblink_id = []
                    for each in speciesdb.dblink_ids_from_gid(gid[0], idtype):
                        dblink_id += each[0]
                    if dblink_id in gene_list:
                        hgene.append(dblink_id)
                        j = gene_list.index(dblink_id)
                        hit_matrix[i,j] = 1
                hit_genes.append(hgene)
                i+=1
            if not 1 in hit_matrix:
                raise exception.GmtdatabaseErr,"Error Message:\nNone of genes in input file(gct file) was matched with genes in gene sets. \nPlease check the gct file(-e), the species(-s), the idtype(-i), the database type(-d)."

        gset_name = np.array(tuple(gset_name))
        gset_des = np.array(gset_des)
        hit_genes = np.array(hit_genes)
        hitsum = hit_matrix.sum(1)
        delindex = np.where((hitsum<min_size) | (hitsum>max_size))[0]
        hit_matrix_filtered = np.delete(hit_matrix,delindex,axis=0)
        hit_genes_filterd = np.delete(hit_genes,delindex,axis=0)
        hit_sum_filtered = np.delete(hitsum, delindex, axis=0)
        gset_name_filtered = np.delete(np.array(gset_name),delindex,axis=0)
        gset_des_filtered = np.delete(gset_des,delindex,axis=0)
        if len(hit_matrix_filtered)==0 and len(gset_name_filtered)==0 and len(gset_des_filtered)==0 and len(hit_genes_filterd)==0 :
            raise exception.GmtvalueErr,"Error Message:\nAll gene sets"+repr(len(hit_matrix))+"have been filtered. \nPlease check the threshold and ceil of gene set size (values of min_size and max_size). "
    except ValueError,e:
        sys.exit(e)
    return gset_num, hit_matrix_filtered, hit_genes_filterd, hit_sum_filtered, gset_name_filtered, gset_des_filtered


def read_cls_file(cls):
    cls_file = open(cls,'r')
    rec = cls_file.readlines()
    try :
        if len(rec)!=3:
            raise exception.ClsformatErr,"Error Message:\nBad phenotype data formats. Please check format of cls file, only \".cls\" format file with a pair of categorical phenotype labels is allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats  "

        line0 = rec[0].rstrip().split()
        if int(line0[1])!=2 :
            raise exception.ClsvalueErr,"Error Message: \nSorry, only two categorical phenotype labels are allowed."
        sample_num = int(line0[0])
        line1 = rec[1].rstrip().split()
        phnameA = line1[1]
        phnameB = line1[2]
        line2 = rec[2].rstrip().split()
        if len(line2)!=sample_num:
            raise exception.ClsvalueErr,"Error Message: \nThe total of phenotype labels doesn't match with the total of samples as you described at the first line. Please check this cls file."
        assert phnameA in line2
#            raise exception.ClsvalueErr,"Error Message:\nThe type name in 3rd line must be same with 2nd line."
        labels = []
        sample0 = 0
        sample1 = 0
        i = 0
        for i in range(sample_num):
            if line2[i]==phnameA:
                labels.append(0)
                sample0 +=1
            else:
                labels.append(1)
                sample1 +=1
#            else :
#                raise exception.ClsvalueErr,"Error Message:\nSorry, we don't support more than 2 phenotype labels."
        cls_file.close()

    except ValueError,e:
        sys.exit(e)
        cls_file.close()
    return labels,phnameA,phnameB,sample_num,sample0,sample1

def read_txt_file(txt, sample_num):
    txt_file=open(txt,'r')
    title = txt_file.readline()
    cols = title.rstrip().split('\t')
    if len(cols)-2 != sample_num :
        raise exception.FormatErr, "Error Message: \nThe total of samples does not match with total of samples in CLS file.\n"
    try:
        lines = txt_file.readlines()
        gene_num = len(lines)
        counter = 0
        expr_data = np.array([[0.0 for i in range(sample_num)] for j in range(gene_num)])
        gene_list = []
        for one in lines:
            words = one.rstrip().split('\t')
            expr_data[counter] = [float(x) for x in words[2:]]
            # if('/' in words[0]):
            #    words[0] = words[0].split(' /// ')[0]
            gene_list.append(words[0])
            counter += 1
        txt_file.close()
#        gene_list = tuple(gene_list)

    except IOError, e:
        sys.exit(e)
    return gene_list, gene_num, expr_data,title


def read_gct_file(gct,sample_num):

    gct_file = open(gct,'r')
    try:
        gct_file.readline()
        num_info = gct_file.readline().rstrip().split('\t')
        if int(num_info[1])!=sample_num:
            raise exception.GctvalueErr,"Error Message: \nThe total samples in gct file doen't match with the total samples in cls file.Please check gct file and cls file."
        gene_num = int(num_info[0])
        title = gct_file.readline()
        header = title.rstrip().split('\t')
        if len(header)!=sample_num+2 and not header[0] in ('NAME','Name','name') :
            raise exception.GctformatErr,"Error Message: \nBad expression data format. Only TXT and GCT format is allowed. \nAbout file format, please click www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats "

        lines = gct_file.readlines()
        if len(lines)!=gene_num :
            raise exception.GctvalueErr,"Error Message: \nThe total gene probes doesn't match with information described in the second line."

        counter = 0
        expr_data = np.array([[0.0 for i in range(sample_num)] for j in range(gene_num)])
        gene_list = []
        for one in lines:
            words = one.rstrip().split('\t')
            expr_data[counter] = [float(x) for x in words[2:]]
            #if('/' in words[0]):
            #    words[0] = words[0].split(' /// ')[0]
            gene_list.append(words[0])
            counter +=1
        gct_file.close()
#        gene_list = tuple(gene_list)
    except IOError,e:
        sys.exit(e)
    return gene_list,gene_num,expr_data,title

def pick_genes(gene_list,gene_num,expr_data, speciesdb, idtype, outdir,title):
    # this function is for picking up genes that are in KOBAS database.
    # use before start compute
    dbGenes = speciesdb.all_genes(idtype)
    rm_idxes = []
    genes_abd = open(outdir+'/Genes_abandoned.txt', 'w')
    genes_anls = open(outdir+'/Genes_Used_In_Analysis.txt','w')
    genes_anls.write(title)
    genes_abd.write(title)
    for idx,gene in enumerate(gene_list):
        if (gene,) in dbGenes:
            gid = speciesdb.gids_from_dblink_id(gene,idtype).fetchone()[0]
            annot = speciesdb.name_from_gid(gid)
            genes_anls.write(gene + '\t' + annot + '\t'+ ('\t').join(map(str,expr_data[idx,])) + '\n')
        else:
            rm_idxes.append(idx)
            genes_abd.write(gene + '\tna\t' + ('\t').join(map(str,expr_data[idx,])) + '\n' ) 
    gene_list = [j for i,j in enumerate(gene_list) if not i in rm_idxes]
    expr_data = np.delete(expr_data, rm_idxes, axis=0)
    if len(gene_list) < gene_num :
        warnings.warn("Warning Messgae: \nSome of your genes can not be annotated, so we removed them for analysis. The genes for analysis are writed in file \" Genes_Used_In_Analysis.txt\" in the output directory. The annotation of each gene according to KEGG GENE is put in 2nd column of this file. The genes abandoned are writed in file  \"Genes_abandoned.txt\". \n") 
    if len(gene_list) != len(expr_data):
        raise Exception, "Error Message: The gene num is not equal to expr_data length\n "
    else:
        gene_num = len(gene_list)
    if len(gene_list) == 0 :
        raise Exception, "Error Message: None of genes are annotated. Please check -e (expression file) -s (species) -i (idtype) \n"
    
    return gene_list, gene_num, expr_data

def rank_pro(lb,md, expr_data,sample0,sample1):


    index_0 = []
    index_1 = []
    for i,j in enumerate(lb):   #abstract index as i, content as j
        if j==0:
            index_0.append(i)

        else :
            index_1.append(i)

    expr_0 = expr_data[:,index_0] # : abstract all rows in column index_0 includes. [row,column]
    expr_1 = expr_data[:,index_1]
    mean_0 = expr_0.mean(1) # axis 1 is meaning get average from column adds
    mean_1 = expr_1.mean(1)
    std_0 = expr_0.std(1)
    std_1 = expr_1.std(1)

    if (md == 'snr'):
        s2n = (mean_0 - mean_1)/(std_0 + std_1)
        sort_gene_index = np.argsort(s2n,0)[::-1].T  #this step get index after sorted, then use this index to get gene list from gene_name
        sort_r = np.sort(s2n,0)[::-1].T #this step get s2n value after sorted

    if (md == 'ttest'):
        a = mean_0 - mean_1
        s0 = np.square(std_0)
        s1 = np.square(std_1)
        b = np.sqrt(s0/sample0 + s1/sample1)
        ttest = a/b
        sort_gene_index = np.argsort(ttest,0)[::-1].T
        sort_r = np.sort(ttest,0)[::-1].T
    return  sort_r, sort_gene_index # NOTE: sort_r is 1*gene_num matrix

def ES_all(sort_r, sort_gene_index, hit_matrix_filtered, weighted_score_type, gene_num):

    hitm = hit_matrix_filtered[:,sort_gene_index]
    missm = hitm - 1
    sort_arr = np.array([sort_r for i in range(len(hitm))])

    if weighted_score_type == 0 :
        tmp = hitm
    if weighted_score_type == 1 :
        tmp = np.absolute(sort_arr)*hitm
    if weighted_score_type == 2 :
        tmp = sort_arr ** 2 * hitm
    else :
        tmp = np.absolute(sort_arr) ** weighted_score_type * hitm

    NR = np.sum(tmp, axis=1)
    hit_score = np.array(map(np.divide,tmp,NR))
    miss_score = 1.0/(gene_num - len(hitm)) *missm
    pre_score = hit_score + miss_score

    RES = np.cumsum(pre_score,axis=1)
    es_idx = lambda x: np.where(abs(x.max()) > abs(x.min()), (x.max(), x.argmax()), (x.min(), x.argmin()))
    re = np.array(map(es_idx,RES))
    es = np.array(re[:,0])
    idx = np.array(re[:,1],dtype = np.int64)
    return es,idx,RES


def ES_for_permutation(lb, md, sample0, sample1, hit_matrix_filtered, weighted_score_type, expr_data, gene_num):

    index_0 = []
    index_1 = []
    for i,j in enumerate(lb):   #abstract index as i, content as j
        if j==0:
            index_0.append(i)
        else :
            index_1.append(i)

    expr_0 = expr_data[:,index_0] # : abstract all rows in column index_0 includes. [row,column]
    expr_1 = expr_data[:,index_1]
    mean_0 = expr_0.mean(1) # axis 1 is meaning get average from column adds
    mean_1 = expr_1.mean(1)
    std_0 = expr_0.std(1)
    std_1 = expr_1.std(1)

    if (md == 'snr'):
        s2n = (mean_0 - mean_1)/(std_0 + std_1)
        sort_gene_index = np.argsort(s2n,0)[::-1].T  #this step get index after sorted, then use this index to get gene list from gene_name
        sort_r = np.sort(s2n,0)[::-1].T #this step get s2n value after sorted

    if (md == 'ttest'):
        a = mean_0 - mean_1
        s0 = np.square(std_0)
        s1 = np.square(std_1)
        b = np.sqrt(s0/sample0 + s1/sample1)
        ttest = a/b
        sort_gene_index = np.argsort(ttest,0)[::-1].T
        sort_r = np.sort(ttest,0)[::-1].T

    hitm = hit_matrix_filtered[:,sort_gene_index]
    missm = hitm - 1
    sort_arr = np.array([sort_r for i in range(len(hitm))])

    if weighted_score_type == 0 :
        tmp = hitm
    if weighted_score_type == 1 :
        tmp = np.absolute(sort_arr)*hitm
    if weighted_score_type == 2 :
        tmp = sort_arr ** 2 * hitm
    else :
        tmp = np.absolute(sort_arr) ** weighted_score_type * hitm

    NR = np.sum(tmp, axis=1)
    hit_score = np.array(map(np.divide,tmp,NR))
    miss_score = 1.0/(gene_num - len(hitm)) *missm
    RES = np.cumsum(hit_score + miss_score,axis=1)
    get_es = lambda x: np.where(abs(x.max()) > abs(x.min()), x.max(), x.min())
    es = np.array(map(get_es,RES))
    return es


def ES_null(lb, times, method, sample0,sample1, hit_matrix_filtered, weighted_score_type,expr_data,gene_num):
    lb_matrix = np.array([lb for i in range(times)])
    ran_labels = np.array(map(lambda x: np.random.permutation(x), lb_matrix))
    def_get_es_null = lambda x: ES_for_permutation(x,method, sample0, sample1, hit_matrix_filtered,weighted_score_type,expr_data,gene_num)
    es_null = np.array(map(def_get_es_null,ran_labels)).T
    return es_null


def nominal_p(es,es_null):

    es_all = np.column_stack((es_null,es))
    def_pval = lambda x: sum(x[:-1] >= x[-1])/float(sum(x[:-1] >=0)) if (x[-1]>=0) else sum(x[:-1] <= x[-1])/float(sum(x[:-1] <=0))
    #def_pval = lambda x: where(x[-1]>=0, sum(x[:-1] >= x[-1])/float(sum(x[:-1]) >=0), sum(x[:-1] <= x[-1])/float(sum(x[:-1]) <=0))
    r = map(def_pval,es_all)
    pval = np.array([r]).T     #m*1 array  m: num of filter gene sets
    return pval

def normalized(es, es_null):

    def_mean_pos = lambda x : np.mean(x[x>=0])
    def_mean_neg = lambda x : np.mean(abs(x[x<=0]))
    def_nor = lambda x : np.where(x[:-2]>=0, x[:-2]/x[-2], x[:-2]/x[-1])

    mean_p = np.array(map(def_mean_pos, es_null))
    mean_p = mean_p.reshape(len(mean_p),1)          # shape=(m,1)
    mean_n = np.array(map(def_mean_neg, es_null))
    mean_n = mean_n.reshape(len(mean_n),1)

    es_null_mean = np.column_stack((es_null,mean_p,mean_n))
    nes_null = np.array(map(def_nor, es_null_mean))

    es_mean = np.column_stack((es,mean_p, mean_n))
    nes_obs = np.array(map(def_nor, es_mean))
    return nes_obs,nes_null

def fdr_cal(nes_obs, nes_null):

    nullmore0 = np.sum(nes_null >= 0)
    nullless0 = np.sum(nes_null <= 0)

    obsmore0 = np.sum(nes_obs >= 0)
    obsless0 = np.sum(nes_obs <= 0)

    def_top = lambda x: np.where(x>=0, np.sum(nes_null >= x)/float(nullmore0) if (nullmore0>0) else 0,
                              np.sum(nes_null <= x)/float(nullless0) if(nullless0>0) else 0)

    def_down = lambda x: np.where(x>=0, np.sum(nes>= x)/float(obsmore0) if(obsmore0>0) else 0,
                               np.sum(nes <= x)/float(obsless0) if(obsless0>0) else 0)

    #def_top = lambda x : where(x>=0, sum(nes_null >= x)/float(nullmore0), sum(nes_null <= x)/float(nullless0))
    #def_down = lambda x : where(x>=0, sum(nes>= x)/float(obsmore0), sum(nes <= x)/float(obsless0))

    top = np.array(map(def_top, nes_obs))
    nes = nes_obs.copy()
    down = np.array(map(def_down, nes_obs))
    cal = np.column_stack((top,down))
    def_fdr = lambda x : np.where(x[0]/x[1] <1, x[0]/x[1], 1)
    fdr = np.array(map(def_fdr,cal))
    return fdr

def output_set(result_path, phnameA, phnameB, gset_name_filtered, gset_des_filtered, hit_matrix_filtered, ES, NES, nompval, FDR):

    gset = {'GeneSets': gset_name_filtered}
    gsetdf = pd.DataFrame(gset)

    gsetdf2 = pd.DataFrame(gset_des_filtered,columns=['Description'])
    # gene set description
    hits_num = hit_matrix_filtered.sum(1)
    hits_numdf = pd.DataFrame(hits_num,columns=["Size"], dtype = 'int32')

    esdf = pd.DataFrame(data=ES, columns=['ES'])
    nesdf = pd.DataFrame(data=NES, columns=['NES'])
    pvaldf = pd.DataFrame(nompval,columns=['Nominal p-value'])
    fdrdf = pd.DataFrame(FDR, columns=['FDR'])

    outset = pd.concat([gsetdf,gsetdf2, hits_numdf, esdf,nesdf,pvaldf,fdrdf], axis = 1)

    type0 = outset[outset['ES']>0].sort_values(by = 'NES', ascending=False)
    type1 = outset[outset['ES']<0].sort_values(by = 'NES', ascending=True)

    path0 = result_path +'/'+ phnameA+ '_Enrichment_Gene_Sets'
    path1 = result_path +'/'+ phnameB+ '_Enrichment_Gene_Sets'

    type0.to_csv(path_or_buf=path0, sep='\t', index=False, mode='a')
    type1.to_csv(path_or_buf=path1, sep='\t', index=False, mode='a')

    return type0,type1

def gene_info(result_path, gene_list, phnameA, phnameB, type0, type1, hit_matrix_filtered, sort_r, sort_gene_index, es_idx, RES ):

    hitm = hit_matrix_filtered[:,sort_gene_index]   # dim: gset_num_flt * gene_num
    glist = np.array(gene_list).T[sort_gene_index]  # dim:  0* gene_num

    if len(type0) <20 :
        type0_num = len(type0)
    else :
        type0_num = 20

    if len(type1) < 20:
        type1_num = len(type1)
    else:
        type1_num = 20



    for idx in type0.index[:type0_num]:
        hit_genes = glist[np.where(hitm[idx,:]==1)]
        core_genes = glist[:es_idx[idx]+1][np.where(hitm[idx,:es_idx[idx]+1]==1)]
        res = RES[idx,np.where(hitm[idx,:]==1)].T
        rank = np.array(np.where(hitm[idx,:]==1)).T+1
        rank_score = sort_r[np.where(hitm[idx,:]==1)]
        core_info = np.hstack((['Yes' for i in range(len(core_genes))],['No' for j in range(len(hit_genes)-len(core_genes))]))
        ph0_out = np.column_stack((hit_genes,core_info, rank, rank_score, res))
        ph0df = pd.DataFrame(ph0_out, columns=['Genes', 'Core Enrichment', 'Rank in Gene List', 'Rank Score', 'RES'])
        path = result_path+ '/' + phnameA+'_'+ type0[type0.index==idx]['GeneSets'].values[0] + '_detail_info'
        ph0df.to_csv(path_or_buf=path,sep='\t', index=False)

    for idx in type1.index[:type1_num]:
        hit_genes = glist[np.where(hitm[idx,:]==1)]
        core_genes = glist[es_idx[idx]:][np.where(hitm[idx,es_idx[idx]:]==1)]
        res = RES[idx,np.where(hitm[idx,:]==1)].T
        rank = np.array(np.where(hitm[idx,:]==1)).T+1
        rank_score = sort_r[np.where(hitm[idx,:]==1)]
        core_info = np.hstack((['No' for j in range(len(hit_genes)-len(core_genes))],['Yes' for i in range(len(core_genes))]))
        ph1_out = np.column_stack((hit_genes,core_info, rank, rank_score, res))
        ph1df = pd.DataFrame(ph1_out, columns=['Genes', 'Core Enrichment', 'Rank in Gene List', 'Rank Score', 'RES'])
        path = result_path+ '/'+ phnameB + '_'+ type1[type1.index==idx]['GeneSets'].values[0] + '_detail_info'
        ph1df.to_csv(path_or_buf=path,sep='\t', index=False)


def plots_RES(result_path, phnameA, phnameB, type0,type1,RES, es_idx, hit_matrix_filtered,sort_gene_index, sort_r):

    ph0pdf = PdfPages(result_path +'/'+phnameA+'.pdf')
    ph1pdf = PdfPages(result_path +'/'+phnameB+'.pdf')

    ph0_setidx = type0.index[:20]
    ph1_setidx = type1.index[:20]
    
    res = pd.DataFrame(RES).T
    hitm = hit_matrix_filtered[:,sort_gene_index]
    
    if len(type0) <20 :
        type0_num = len(type0)
    else :
        type0_num = 20    

    if len(type1) < 20:
        type1_num = len(type1)
    else:
        type1_num = 20 
    
    for i in range(type0_num):
        fig = plt.figure(figsize=(8,6))
        plt.axes([0.025,0.025,0.95,0.95])
        gs = gridspec.GridSpec(3,1,height_ratios=[3,1,2])
        plt.subplots_adjust(wspace=0,hspace=0)

        ax1 = plt.subplot(gs[0])
        data = res[ph0_setidx[i]]
        ymin,ymax = data.min(), data.max()
        dy = (ymax - ymin)*0.2
        ax1.set_ylim(ymin-dy,ymax+dy)
        ax1.set_xlim(0,int(1.05*len(data)))
        ax1.grid(which='major',axis='x', linewidth = 0.75, linestyle='-', color = '0.75')
        ax1.grid(which='major',axis='y', linewidth = 0.75, linestyle='-', color = '0.75')
        ax1.plot(data,color='red',label='Running Enrichment Score')
        t = es_idx[type0.index[i]]
        ax1.annotate("ES score", xy=(t,data[t]), xytext = (t+500,1.1*data[t]),
                    arrowprops=dict(facecolor='black', shrink=0.1, width=0.8,headwidth=3))
        ax1.set_ylabel('Running Enrichment Score (RES)')
        plt.legend(loc='upper right', frameon=False)
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = plt.subplot(gs[1],sharex=ax1)
        plt.ylim(0,1)
        ax2.set_yticks([])
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.plot(hitm[ph0_setidx[i]],color='green', label='Hits')

        ax3 = plt.subplot(gs[2],sharex=ax1)
        rmin , rmax = sort_r.min(), sort_r.max()
        dr = 0.1*(rmax-rmin)
        plt.ylim(rmin-dr,rmax+dr)
        ax3.grid(which='major',axis='x', linewidth = 0.75, linestyle='-', color = '0.75')
        ax3.grid(which='major',axis='y', linewidth = 0.75, linestyle='-', color = '0.75')
        ax3.plot(sort_r,label='Ranking metric scores')
        plt.fill_between(data.index,sort_r,0,alpha=.5,color='grey')
        ax3.set_ylabel('Ranking correlationship value')
        ax3.set_xlabel('Rank in Ordered Dataset')
        plt.legend(loc='upper right', frameon=False)
        fig.suptitle(phnameA + ' ' +type0.iloc[i]['GeneSets']+ ' ('+type0.iloc[i]['Description'] +')'+ ' GSEA_Result',size=15)
        ph0pdf.savefig()
        plt.close()
    ph0pdf.close()

    for i in range(type1_num):
        fig = plt.figure(figsize=(8,6))
        plt.axes([0.025,0.025,0.95,0.95])
        gs = gridspec.GridSpec(3,1,height_ratios=[3,1,2])
        plt.subplots_adjust(wspace=0,hspace=0)

        ax1 = plt.subplot(gs[0])
        data = res[ph1_setidx[i]]
        ymin,ymax = data.min(), data.max()
        dy = (ymax - ymin)*0.2
        ax1.set_ylim(ymin-dy,ymax+dy)
        ax1.set_xlim(0,int(1.05*len(data)))
        ax1.grid(which='major',axis='x', linewidth = 0.75, linestyle='-', color = '0.75')
        ax1.grid(which='major',axis='y', linewidth = 0.75, linestyle='-', color = '0.75')
        ax1.plot(data,color='red',label='Running Enrichment Score')
        t = es_idx[type1.index[i]]
        ax1.annotate("ES score", xy=(t,data[t]), xytext = (t+500,1.1*data[t]),
                    arrowprops=dict(facecolor='black', shrink=0.1, width=0.8,headwidth=3))
        ax1.set_ylabel('Running Enrichment Score (RES)')
        plt.legend(loc='upper right', frameon=False)
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2 = plt.subplot(gs[1],sharex=ax1)
        plt.ylim(0,1)
        ax2.set_yticks([])
        plt.setp(ax2.get_xticklabels(), visible=False)
        ax2.plot(hitm[ph1_setidx[i]],color='green', label='Hits')

        ax3 = plt.subplot(gs[2],sharex=ax1)
        rmin , rmax = sort_r.min(), sort_r.max()
        dr = 0.1*(rmax-rmin)
        plt.ylim(rmin-dr,rmax+dr)
        ax3.grid(which='major',axis='x', linewidth = 0.75, linestyle='-', color = '0.75')
        ax3.grid(which='major',axis='y', linewidth = 0.75, linestyle='-', color = '0.75')
        ax3.plot(sort_r,label='Ranking metric scores')
        plt.fill_between(data.index,sort_r,0,alpha=.5,color='grey')
        ax3.set_ylabel('Ranking correlationship value')
        ax3.set_xlabel('Rank in Ordered Dataset')
        plt.legend(loc='upper right', frameon=False)
        fig.suptitle( phnameB + ' ' +type1.iloc[i]['GeneSets']+ ' ('+type1.iloc[i]['Description'] +')'+' GSEA_Result',size=15)
        ph1pdf.savefig()
        plt.close()
    ph1pdf.close()

