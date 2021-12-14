#!/usr/bin/env python

def list_to_tuple(rows):
    nrows = set()
    for row in rows:
        nrows.add(tuple(row))

    return list(nrows)

def dict_to_list(rows):
    nrows = set()
    for row in rows.items():
        nrows.add(tuple(row[1] + [row[0]]))

    return list(nrows)

def dict_to_list_reorder(rows):
    nrows = set()
    for row in rows.items():
        row[1].insert(2, row[0])
        nrows.add(tuple(row[1]))

    return list(nrows)

def split_values(rows):
    nrows = set()
    for row in rows.items():
        for value in row[1]:
            nrows.add((row[0], value))

    return list(nrows)

def combine_dicts(frows, srows):
    nrows = set()
    for frow in frows.items():
        for fvalue in frow[1]:
            if srows.has_key(fvalue):
                for svalue in srows[fvalue]:
                    nrows.add((frow[0], svalue))

    return list(nrows)
