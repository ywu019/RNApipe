from StringIO import StringIO

from Bio.Blast import NCBIXML, Applications

import sys

from kobas import config, dbutils

class Annotation(object):

    def __init__(self, query, links):
        self.query = query
        self.links = links

    def has_links(self):
        return len(self.links) != 0

def parse_annot(line):
    links = []
    try:
        query, slinks = line[:-1].split('\t')
    except ValueError:
	sys.exit('Error message:\nBad format of the input file. You may do annotation first or check your input file. The input file should be the output of annotate.')
    if slinks != 'None':
        slinks = slinks.split('!')
        for slink in slinks:
            links.append(slink.split('|')[0])
    return Annotation(query, tuple(links))

def is_info(line):
    return (line.strip() != '' and line[0] != '#')

def get_args(handle):
    first_line = handle.readline()
    abbr = first_line.split()[0][2:]
    second_line = handle.readline()
    if second_line.startswith('##Method: ID mapping'):
        idmapping = second_line.split()[-1][1:-1]
    else:
        idmapping = 'entrez_gene_id'
    third_line = handle.readline()
    if third_line.startswith('##Summary'):
        fgnum = third_line.split('\t')[1].split()[0]
    return abbr, idmapping, fgnum


class Iterator(object):

    def __init__(self, handle):
        self._handle = handle
        self.flag = 1

    def __iter__(self):
        return self

    def next(self):
        while True:
            line = self._handle.next()
            if line.startswith('-' * 20):
                self.flag = 0
            if (self.flag == 1 and is_info(line)):
                return parse_annot(line)

class Annotator(object):

    def __init__(self, reader, selector):
        self.reader, self.selector = reader, selector

    def annotate(self, verbose=True):
        for record in self.reader:
            if not record:
                break
            yield self.selector.select(record)

class Reader(object):

    def __iter__(self):
        raise NotImplementedError

class BlastoutXMLReader(Reader):

    def __init__(self, handle):
        self.handle = handle

    def __iter__(self):
        return NCBIXML.parse(self.handle)

class BlastProgReader(BlastoutXMLReader):

    def __init__(self, program, blastcmd, infile, database, nCPUs):
        cline = getattr(Applications, 'Ncbi' + program + 'Commandline')(cmd = blastcmd, query = infile, db = database, outfmt = 5, num_threads = nCPUs)
        blastout, blasterr = cline()
        BlastoutXMLReader.__init__(self, StringIO(blastout))

class BlastoutTabReader(Reader):

    def __init__(self, handle):
        records = {}
        for line in handle:
            if (line[0] != '#' and line.split() != []):
                record = line.split('\t')
                records.setdefault(record[0], []).append(record[1:])
        self.parse = records.items().__iter__()

    def __iter__(self):
        return self.parse

class IdMappingReader(file):

    def __init__(self, speciesdb, species, infile, dbtype):
        self.speciesdb = speciesdb
        self.species = species
        file.__init__(self, infile)
        self.dbtype = dbtype

    def next(self):
        dblink_id = file.next(self).strip()
        while dblink_id == '':    ##remove blank lines
            dblink_id = file.next(self).strip()
        if self.species == 'ko':
            return (dblink_id, self.speciesdb.koids_from_dblink_id(dblink_id, self.dbtype))
        else:
            return (dblink_id, self.speciesdb.gids_from_dblink_id(dblink_id, self.dbtype))

class Selector(object):

    def select(self, record):
        raise NotImplementedError

class BlastoutXMLSelector(Selector):

    def __init__(self, speciesdb, species, cutoffs):
        self.speciesdb, self.species, self.cutoffs = speciesdb, species, cutoffs

    def select(self, record):
        ret, rank = Annotation(record.query.split()[0], set()), 1
        for alignment in record.alignments:
            if rank > self.cutoffs[0] or alignment.hsps[0].expect > self.cutoffs[1] or \
                (float(alignment.hsps[0].align_length) / alignment.length) < self.cutoffs[2]:
                break
            gid = str(alignment.hit_def)
            if self.species == 'ko':
                koids = self.speciesdb.koids_from_gid(gid)
                for koid in koids:
                    name = self.speciesdb.name_from_koid(koid[0])
                    ret.links.add((koid[0], name))
                    return ret
            else:
                name = self.speciesdb.name_from_gid(gid)
                ret.links.add((gid, name))
                return ret    ##rank is only used for KO
            rank += 1
        return ret

class BlastoutTabSelector(Selector):

    def __init__(self, speciesdb, species, cutoffs):
        self.speciesdb, self.species, self.cutoffs = speciesdb, species, cutoffs

    def select(self, record):
        ret, rank = Annotation(record[0], set()), 1
        for alignment in record[1]:
            if rank > self.cutoffs[0] or float(alignment[9]) > self.cutoffs[1]:
                break
            gid = alignment[0]
            if self.species == 'ko':
                koids = self.speciesdb.koids_from_gid(gid)
                for koid in koids:
                    name = self.speciesdb.name_from_koid(koid[0])
                    ret.links.add((koid[0], name))
                    return ret
            else:
                name = self.speciesdb.name_from_gid(gid)
                ret.links.add((gid, name))
                return ret    ##rank is only used for KO
            rank += 1
        return ret

class IdMappingSelector(Selector):

    def __init__(self, speciesdb, species):
        self.speciesdb, self.species, = speciesdb, species

    def select(self, record):
        ret = Annotation(record[0], set())
        if self.species == 'ko':
            for koid in record[1]:
                name = self.speciesdb.name_from_koid(koid[0])
                ret.links.add((koid[0], name))
        else:
            for gid in record[1]:
                name = self.speciesdb.name_from_gid(gid[0])
                ret.links.add((gid[0], name))
        return ret
