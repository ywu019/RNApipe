#!/usr/bin/env python

def parse(handle):
    symbol_geneids, synonym_geneids = {}, {}

##tax_id\tGeneID\tSymbol\tLocusTag\tSynonyms\t
    for line in handle:
        geneid, symbol, tmp, synonyms = line.split('\t')[1:5]

        symbol_geneids.setdefault(symbol, set()).add(geneid)

        if synonyms != '-':
            synonyms = synonyms.split('|')
            for synonym in synonyms:
                synonym_geneids.setdefault(synonym, set()).add(geneid)

    for symbol_geneid in symbol_geneids.items():
        if len(symbol_geneid[1]) > 1:
            symbol_geneids.pop(symbol_geneid[0])
        else:
            symbol_geneids[symbol_geneid[0]] = list(symbol_geneid[1])[0]

    for synonym_geneid in synonym_geneids.items():
        if len(synonym_geneid[1]) > 1:
            synonym_geneids.pop(synonym_geneid[0])
        else:
            synonym_geneids[synonym_geneid[0]] = list(synonym_geneid[1])[0]

    synonym_geneids.update(symbol_geneids)

    return synonym_geneids

if __name__ == '__main__':
    import sys
    from pprint import pprint

    symbol_geneids = parse(open(sys.argv[1]))
    pprint(symbol_geneids)
