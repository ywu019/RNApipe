#!/usr/bin/python

import os, sys, re, tempfile
from optparse import OptionParser

from kobas import annot, config, dbutils, exception, fasta, output, discover
def config_option():
    '''configure options and arguments for annotate.py
    '''
    usage = 'Usage: run_kobas.py [-l] -i infile [-t intype] -s species [-E evalue] [-R rank] [-N nCPUs] [-C coverage] [-Z ortholog] [-b bgfile] [-d database] [-m method] [-n fdr] [-o outfile] [-c cutoff] [-k kobas_home] [-v blast_home] [-y blastdb] [-q kobasdb] [-p blastp] [-x blastx]'
    p = OptionParser(usage)
    # options for information
    p.add_option(
        '-l', '--list', dest = 'list', action = 'store_true',
        help = 'list available species, or list available databases for specific species')
    # basic options
    p.add_option(
        '-i', '--infile', dest = 'infile', action = 'store',
        type = 'string', help = 'input data file')
    p.add_option(
        '-t', '--intype', dest = 'intype', default = 'fasta:pro', action = 'store',
        type = 'string', help = 'input type (%s, blastout:xml, blastout:tab, %s), default fasta:pro' %
            (', '.join(PROGRAMS.keys()), ', '.join(DBLINKS.keys())))
    p.add_option(
        '-s', '--species', dest = 'species', action = 'store',
        type = 'string', help = 'species abbreviation (for example: ko for KEGG Orthology, hsa for Homo sapiens, mmu for Mus musculus, dme for Drosophila melanogaster, ath for Arabidopsis thaliana, sce for Saccharomyces cerevisiae and eco for Escherichia coli K-12 MG1655)')
    # options for blast and parsing blast result
    p.add_option(
        '-E', '--evalue', dest = 'evalue', default = 1e-5, action = 'store',
        type = 'float', help = 'expect threshold for BLAST, default 1e-5')
    p.add_option(
        '-R', '--rank', dest = 'rank', default = 5, action = 'store',
        type = 'int', help = 'rank cutoff for valid hits from BLAST result, default 5')
    p.add_option(
        '-N', '--nCPUs', dest = 'nCPUs', default = 1, action = 'store',
        type = 'int', help = 'number of CPUs to be used by BLAST, default 1')
    # option for reviewer
    p.add_option(
        '-C', '--coverage', dest = 'coverage', default = 0.0, action = 'store',
        type = 'float', help = 'subject coverage cutoff for BLAST, default 0')
    p.add_option(
        '-Z', '--ortholog', dest = 'ortholog', default = 'NO', action = 'store',
        type = 'string', help = 'whether only use ortholog for cross-species annotation or not, default NO (If only use ortholog, give species abbr)')
    p.add_option(
        '-b', '--bgfile', dest = 'bgfile', default = 'same', action = 'store',
        type = 'string', help = 'background file, the output of annotate (3 or 4-letter file name is not allowed), or species abbreviation (for example: hsa for Homo sapiens, mmu for Mus musculus, dme for Drosophila melanogaster, ath for Arabidopsis thaliana, sce for Saccharomyces cerevisiae and eco for Escherichia coli K-12 MG1655), default same species as annotate')
    p.add_option(
        '-d', '--db', dest = 'db', default = 'K/R/B/p/o/k/N/G/S', action = 'store',
        type = 'string', help = 'databases for selection, 1-letter abbreviation separated by "/": K for KEGG PATHWAY, R for Reactome, B for BioCyc, p for PANTHER, o for OMIM, k for KEGG DISEASE, N for NHGRI GWAS Catalog and G for Gene Ontology, S for Gene Ontology Slim, default K/R/B/p/o/k/N/G/S')
    p.add_option(
        '-m', '--method', dest = 'method', default = 'h', action = 'store',
        type = 'string', help = "choose statistical test method: b for binomial test, c for chi-square test, h for hypergeometric test / Fisher's exact test, and x for frequency list, default hypergeometric test / Fisher's exact test")
    p.add_option(
        '-n', '--fdr', dest = 'fdr', default = 'BH', action = 'store',
        type = 'string', help = 'choose false discovery rate (FDR) correction method: BH for Benjamini and Hochberg, BY for Benjamini and Yekutieli, QVALUE, and None, default BH')
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for identification result, default stdout')
    p.add_option(
        '-c', '--cutoff', dest = 'cutoff', default = 5, action = 'store',
        type = 'int', help = 'the gene number in a term is not less than the cutoff, default 5')
## add to replace kobasrc
    p.add_option(
        '-k', '--kobashome', dest = 'kobas_home', default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to kobas_home, which is parent directory of sqlite3/ and seq_pep/ , default value is read from ~/.kobasrcwhere you set before running kobas. If you set this parameter, it means you set "kobasdb" and "blastdb" in this following directory. e.g. "-k /home/user/kobas/", means that you set kobasdb = /home/user/kobas/sqlite3/ and blastdb = /home/user/kobas/seq_pep/'
    )
    p.add_option(
        '-v', '--blasthome', dest = 'blast_home', default = '', action = 'store',
        type = 'string',help = 'Optional parameter. To set parent directory of blastx and blastp. If you set this parameter, it means you set "blastx" and "blastp" in this following directory. Default value is read from ~/.kobasrc where you set before running kobas'
    )
    p.add_option(
        '-y','--blastdb', dest = 'blastdb', default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to sep_pep/, default value is read from ~/.kobasrc where you set before running kobas'
    )
    p.add_option(
        '-q','--kobasdb',dest = 'kobasdb',default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to sqlite3/, default value is read from ~/.kobasrc where you set before running kobas, e.g. "-q /kobas_home/sqlite3/"'
    )
    p.add_option(
        '-p', '--blastp', dest = 'blastp', default = '', action = 'store',
        type = 'string', help = 'Optional parameter. To set path to blastp program, default value is read from ~/.kobasrc where you set before running kobas')
    p.add_option(
        '-x','--blastx', dest = 'blastx', default = '', action = 'store',
        type = 'string', help= 'Optional parameter. To set path to  blasx program, default value is read from ~/.kobasrc where you set before running kobas'
    )
    opt, args = p.parse_args()
    return (p, opt, args)

PROGRAMS = {'fasta:pro': 'blastp', 'fasta:nuc': 'blastx'}
DBLINKS = {'id:ncbigene': 'entrez_gene_id', 'id:refseqpro': 'refseq_protein_id', 'id:uniprot': 'uniprotkb_ac', 'id:ensembl': 'ensembl_gene_id'}

METHODS = {'b': 'BinomTest', 'c': 'ChiTest', 'h': 'FisherTest', 'x': 'FreqList'}
FULLNAMES = {'b': 'binomial test', 'c': 'chi-square test', 'h': "hypergeometric test / Fisher's exact test",
             'BH': 'Benjamini and Hochberg', 'BY': 'Benjamini and Yekutieli', 'QVALUE': 'QVALUE'}

opt_parser, opt, args = config_option()

if len(sys.argv) == 1:
    opt_parser.print_help()
    sys.exit(1)

pattern = re.compile('run_kobas')
cline_a = 'python '+ pattern.sub('annotate',sys.argv[0])
cline_i = 'python '+ pattern.sub('identify',sys.argv[0])

keys = lambda y,x:(y+str(x)) if(x) else('')

clinekeys_a = ((' -l') if (opt.list) else (''))+\
              keys(' -i ', opt.infile)+\
              keys(' -s ', opt.species)+\
              keys(' -t ', opt.intype)+\
              keys(' -e ', opt.evalue)+\
              keys(' -r ', opt.rank)+\
              keys(' -n ', opt.nCPUs)+\
              keys(' -c ', opt.coverage)+\
              keys(' -z ', opt.ortholog)+\
              keys(' -k ', opt.kobas_home)+\
              keys(' -v ', opt.blast_home)+\
              keys(' -y ', opt.blastdb)+\
              keys(' -q ', opt.kobasdb)+\
              keys(' -p ', opt.blastp)+\
              keys(' -x ', opt.blastx)

clinekeys_i = keys(' -b ', opt.bgfile)+\
              keys(' -d ', opt.db)+\
              keys(' -m ', opt.method)+\
              keys(' -n ', opt.fdr)+\
              keys(' -o ', opt.outfile)+\
              keys(' -c ', opt.cutoff)+\
              keys(' -k ', opt.kobas_home)+\
              keys(' -v ', opt.blast_home)+\
              keys(' -y ', opt.blastdb)+\
              keys(' -q ', opt.kobasdb)+\
              keys(' -p ', opt.blastp)+\
              keys(' -x ', opt.blastx)

cline_a+=clinekeys_a
cline_i+=clinekeys_i

if(opt.list):
    print os.popen(cline_a).read(),
    sys.exit()

with tempfile.NamedTemporaryFile('wt') as temp1:
    temp1.write(os.popen(cline_a).read())
    temp1.seek(0)
    cline_i+=(' -f '+temp1.name)
    os.system(cline_i)
