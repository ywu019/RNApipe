from Bio import SeqIO

import sys,os

from Bio.Blast import NCBIXML

from exception import *

MaxSeqNo = 500

def verify(handle):
    fasta = SeqIO.parse(handle,'fasta')
#    try:
#        records=list(fasta)
#    	if not records:
#            raise FastaErr, 'Error message:\nBad fasta format. Please check your input file or change the input type.'
#    except FastaErr, e:
#	sys.exit(e)
    try:
        count = 0
        for rec in fasta:
            ##capture the error of only one title without sequence at the end of fasta file, ignoring blank lines at the end.
            if rec.seq.tostring() == '' or rec.id == '':
                raise IOError
            count += 1
            if count > MaxSeqNo:
                pass
                #raise FastaIOError, '''The sequences you have submitted exceed the MAXIMUM LIMIT of %d sequences. Too many sequences may overburden KOBAS 2.0 web server. It is recommended to run it locally using <a href="http://kobas.cbi.pku.edu.cn/download.do">KOBAS 2.0 Standalone programs</a>.''' % MaxSeqNo
        if count == 0:
            raise FastaErr, 'Error message:\nBad fasta format. Please check your input file or change the input type.'
    except FastaErr, e:
        sys.exit(e)
    except FastaIOError, e:
        raise
    except Exception, args:
        seq = locals().get('rec', 'the first')
        if seq.id:
            label = locals().get('rec', 'the first').id
            msg = 'Bad fasta format at " the FIRST " sequence'
        else:
            label = '\n%s\n' % str(seq)
            msg = 'Bad fasta format after " %s " sequence' % label
        raise IOError, msg

def blastoutxml(handle):
    try:
        for record in NCBIXML.parse(handle):
	    test=record
    except ValueError:
	sys.exit('Error message:\nBad blastout xml format. Please check your input file.')

def blastouttab(handle):
    for line in handle:
        if (line[0] != '#' and line.split() != []):
	    lines=line.split('\t')
	    colums=len(lines)
	    try:
	        if not colums==12:
		    raise TabformatErr,'Error message:\nBad blastout tabular format. Please check your input file.'
	    except TabformatErr, e:
                sys.exit(e)
