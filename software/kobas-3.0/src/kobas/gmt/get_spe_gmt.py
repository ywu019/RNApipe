import sys
from kobas import config
kobasrc = config.getrc()

spe = sys.argv[1]
outpath = sys.argv[2]

if len(sys.argv)!= 3:
    print "Usage: python "+ sys.argv[0] +'species outpath \n'
    print "Output: bash scripts for this species of all gmts" 

db_tp = {'K':'P', 'R':'P', 'B':'P', 'p':'P','k':'D','G':'G','S':'G'}

dbs = db_tp.keys()

idtypes = ['ensembl','ncbigene']

cmd = 'python /gpfs/user/aic/kobas-3.0.0/src/kobas/gmt/getGmt.py '+ spe +' '

with open(outpath + '/'+spe + '_gmt_files.sh', 'w') as f :
    for idtype in idtypes:
        for db in dbs:
            cmd2 = cmd+ db_tp[db]+ ':' + db + ' ' + idtype + ' '+kobasrc['gmt'] 
            f.write(cmd2+'\n')

