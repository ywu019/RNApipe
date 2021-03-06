#!/usr/bin/env python

import sys#, time

from optparse import OptionParser

from kobas import annot, config, dbutils, fasta, exception, output

def config_option():
    usage = 'Usage: %prog [-l] -i infile [-t intype] -s species [-o outfile] [-e evalue] [-r rank] [-n nCPUs] [-c coverage] [-z ortholog] [-k kobas_home] [-v blast_home] [-y blastdb] [-q kobasdb] [-p blastp] [-x blastx]'
    p = OptionParser(usage)
    ##options for information
    p.add_option(
        '-l', '--list', dest = 'list', action = 'store_true',
        help = 'list available species, or list available databases for a specific species')
    ##basic options
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
    p.add_option(
        '-o', '--outfile', dest = 'outfile', action = 'store',
        type = 'string', help = 'output file for annotation result, default stdout')
    ##options for blast and parsing blast result
    p.add_option(
        '-e', '--evalue', dest = 'evalue', default = 1e-5, action = 'store',
        type = 'float', help = 'expect threshold for BLAST, default 1e-5')
    p.add_option(
        '-r', '--rank', dest = 'rank', default = 5, action = 'store',
        type = 'int', help = 'rank cutoff for valid hits from BLAST result, default 5')
    p.add_option(
        '-n', '--nCPUs', dest = 'nCPUs', default = 1, action = 'store',
        type = 'int', help = 'number of CPUs to be used by BLAST, default 1')
    ##option for reviewers
    p.add_option(
        '-c', '--coverage', dest = 'coverage', default = 0.0, action = 'store',
        type = 'float', help = 'subject coverage cutoff for BLAST, default 0')
    p.add_option(
        '-z', '--ortholog', dest = 'ortholog', default = 'NO', action = 'store',
        type = 'string', help = 'whether only use orthologs for cross-species annotation or not, default NO (if only use orthologs, please provide the species abbreviation of your input)')
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

opt_parser, opt, args = config_option()

##KOBAS environment configuration
kobasrc = config.getrc()
if opt.kobas_home:
    kobasrc['kobas_home'] = opt.kobas_home
    kobasrc['kobasdb'] = opt.kobas_home + '/sqlite3/'
    kobasrc['blastdb'] = opt.kobas_home + '/seq_pep/'
if opt.blast_home:
    kobasrc['blast_home'] = opt.blast_home
    kobasrc['blastp'] = opt.blast_home + '/blastp/'
    kobasrc['blastx'] = opt.blast_home + '/blastx/'
if opt.blastdb:
    kobasrc['blastdb'] = opt.blastdb +'/'
if opt.kobasdb:
    kobasrc['kobasdb'] = opt.kobasdb +'/'
if opt.blastp:
    kobasrc['blastp'] = opt.blastp +'/'
if opt.blastx:
    kobasrc['blastx'] = opt.blastx +'/'

##open KOBASDB
organismdb = dbutils.KOBASDB(kobasrc['kobasdb'] + 'organism.db')
if opt.species:
    speciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + opt.species + '.db')

if opt.list:
    if opt.species:
        print 'Available databases for %s:' % opt.species
        databases = speciesdb.databases_from_abbr(opt.species)
        for database in databases:
            print '\t'.join(database)
    else:
        print 'Available species: \nko\tKEGG Orthology'
        species = organismdb.organisms(name = True)
        for specie in species:
            print '\t'.join(specie)
    sys.exit()

if opt.infile:
    args.insert(0, opt.infile)
else:
    opt_parser.error('Option -i must be assigned.\n')

if not opt.species:
    opt_parser.error('Option -s must be assigned.\n')

if opt.outfile:
    global old_stdout
    old_stdout = sys.stdout
    sys.stdout = open(opt.outfile, 'w')

if opt.intype in PROGRAMS.keys():
    ##verify fasta file
    try:
        f = open(args[0])
        try:
            fasta.verify(f)
        finally:
            f.close()
    except exception.FastaIOError, msg:
        exception.error(msg)
        sys.exit(1)
    ##key step
    program = PROGRAMS[opt.intype]
    annotator = annot.Annotator(
    reader = annot.BlastProgReader(
        program, kobasrc[program], args[0], kobasrc['blastdb'] + opt.species + '.pep.fasta', opt.nCPUs),
    selector = annot.BlastoutXMLSelector(speciesdb, opt.species, [opt.rank, opt.evalue, opt.coverage]))

elif opt.intype == 'blastout:xml':
    f=open(args[0])
    try:
        fasta.blastoutxml(f)
    finally:
        f.close()
    annotator = annot.Annotator(
        reader = annot.BlastoutXMLReader(open(args[0])),
        selector = annot.BlastoutXMLSelector(speciesdb, opt.species, [opt.rank, opt.evalue, opt.coverage]))
elif opt.intype == 'blastout:tab':
    f=open(args[0])
    try:
        fasta.blastouttab(f)
    finally:
        f.close()
    annotator = annot.Annotator(
        reader = annot.BlastoutTabReader(open(args[0])),
        selector = annot.BlastoutTabSelector(speciesdb, opt.species, [opt.rank, opt.evalue]))
elif opt.intype in DBLINKS.keys():
    dbtype = DBLINKS[opt.intype]
    annotator = annot.Annotator(
        reader = annot.IdMappingReader(speciesdb, opt.species, args[0], dbtype),
        selector = annot.IdMappingSelector(speciesdb, opt.species))
else:
    sys.exit('%s input is not supported yet. Only %s, blastout:xml, blastout:tab, %s are supported.' %
        (opt.intype, ', '.join(PROGRAMS.keys()), ', '.join(DBLINKS.keys())))

##if query ids are same, only keep the first one
items, item_querys = [item for item in annotator.annotate()], []
for item in tuple(items):
    if item.query not in item_querys:
        item_querys.append(item.query)
    else:
        items.remove(item)

##orthology
if (opt.intype in PROGRAMS.keys() and opt.species != 'ko' and opt.ortholog != 'NO'):
    ospeciesdb = dbutils.KOBASDB(kobasrc['kobasdb'] + opt.ortholog + '.db')
    ##run BLAST for the species itself
    annotator = annot.Annotator(
    reader = annot.BlastProgReader(
        program, kobasrc[program], args[0], kobasrc['blastdb'] + opt.ortholog + '.pep.fasta', opt.nCPUs),
    selector = annot.BlastoutXMLSelector(ospeciesdb, opt.ortholog, [opt.rank, opt.evalue, opt.coverage]))
    ##filter with ortholog
    orthologs = {}
    for item in annotator.annotate():
        if not orthologs.has_key(item.query):
            orthologs[item.query] = item.links
    for item in items:
        if not (item.links != set() and orthologs[item.query] != set() and speciesdb.is_ortholog(list(item.links)[0][0], list(orthologs[item.query])[0][0])):
            items[items.index(item)].links = set()

##report annotation result
num_genes_has_annot, num_genes_hasnot_annot = 0, 0

for item in items:
    if item.has_links():
        num_genes_has_annot += 1
    else:
        num_genes_hasnot_annot += 1

##check the number of genes which have been annotated
try:
    if num_genes_has_annot == 0:
        raise
except Exception:
    sys.exit('None of the entries was annotated successfully. Please check the input file(-i), type(-t) and species(-s).')

#gmttime = time.asctime(time.localtime(time.mktime(time.gmtime())-3600*8))
#print '##Time\t%s (UTC)' % gmttime
#print '##Program version: December 1st, 2013  Database version: March 20th, 2014'
if opt.species == 'ko':
    print '##ko\tKEGG Orthology'
else:
    species_name = organismdb.name_from_abbr(opt.species)
    print '##%s\t%s' % (opt.species, species_name)
if opt.intype in ('fasta:pro', 'fasta:nuc', 'blastout:xml', 'blastout:tab'):
    if opt.species == 'ko':
        print '##Method: BLAST\tOptions: evalue <= %s; rank <= %s' % (str(opt.evalue), str(opt.rank))
    else:
        print '##Method: BLAST\tOptions: evalue <= %s' % str(opt.evalue)
elif opt.intype in DBLINKS.keys():
    print '##Method: ID mapping (%s)' % dbtype
print '##Summary:\t%d succeed, %d fail' % (num_genes_has_annot, num_genes_hasnot_annot)
if opt.species == 'ko':
    print '\n#Query\tKO ID|KO name|Hyperlink'
else:
    print '\n#Query\tGene ID|Gene name|Hyperlink'

output.annotate_table(items, opt.species)
output.annotate_text(items, opt.species, speciesdb)
