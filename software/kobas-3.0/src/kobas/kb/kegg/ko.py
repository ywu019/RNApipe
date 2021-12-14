#!/usr/bin/env python

import re

from kobas import utils

def parse(handle):
    kos, pathways, ko_pathways = [], {}, []

    is_pathway = 0
    for line in handle:
        if line.startswith('ENTRY'):
            koid = line.split()[1]
            kos.append([koid, ''])
        if line.startswith('NAME'):
            kos[-1][1] = line[:-1].split(None, 1)[1]

        if line.startswith('PATHWAY'):
            is_pathway = 1
        elif not line.startswith(' '):
            is_pathway = 0

        if is_pathway:
            line = re.sub('PATHWAY', '', line)
            pid, pname = line[:-1].split(None, 1)
            if not pathways.has_key(pid):
                pathways[pid] = pname
            ko_pathways.append((koid, pid))

    kos = utils.list_to_tuple(kos)
    pathways = pathways.items()
    ko_pathways = utils.list_to_tuple(ko_pathways)

    return (kos, pathways, ko_pathways)

if __name__ == '__main__':
    import sys
    from pprint import pprint

    kos, pathways, ko_pathways = parse(open(sys.argv[1]))

    print len(kos), len(pathways), len(ko_pathways)

    k, p = set(), set()
    for kp in ko_pathways:
        k.add(kp[0])
        p.add(kp[1])

    print len(k), len(p)

    pprint(kos[:5])
    pprint(pathways[:5])
    pprint(ko_pathways[:5])
