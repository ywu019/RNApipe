#!/usr/bin/env python

import os, sys
from ConfigParser import ConfigParser

def getrc():
    rc = KOBASRC()

    if os.environ.has_key('HOME'):
        rcfg_file = os.environ['HOME'] + '/.kobasrc'
    else:
        rcfg_file = ''

    rcbg_file = '/etc/kobasrc'

    rcfg_filep = os.access(rcfg_file, os.F_OK)
    rcbg_filep = os.access(rcbg_file, os.F_OK)

    if rcfg_filep and rcbg_filep:
        return rc.merge(rcfg_file, rcbg_file)
    elif rcfg_filep:
        return rc.read(rcfg_file)
    elif rcbg_filep:
        return rc.read(rcbg_file)
    else:
        print 'Error: configuration file does not exist.'
        sys.exit(1)

class KOBASRC:
    def read(self, rc_file):
        rc = {}
        rcparser = ConfigParser()
        rcparser.read(rc_file)

        for option, value in rcparser.defaults().items():
            rc[option] = value

        for section in rcparser.sections():
            for option in rcparser.options(section):
                rc[option] = rcparser.get(section, option)

        return rc

    def merge(self, rcfg_file, rcbg_file):
        '''Options in rcfg_file override the ones in rcbg_file.
        '''
        rcfg = self.read(rcfg_file)
        rcbg = self.read(rcbg_file)
        rcbg.update(rcfg)
        return rcbg

if __name__ == '__main__':
    from pprint import pprint

    pprint(getrc())
