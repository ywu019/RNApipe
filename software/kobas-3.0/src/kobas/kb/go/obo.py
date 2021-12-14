#!/usr/bin/env python

class GoTerm:

    def __init__(self):
        self.relations = set()
        self.obsolete = False

    def set_goid(self, goid):
        self.goid = goid

    def set_name(self, name):
        self.name = name

    def add_relation(self, relation_type, goid):
        self.relations.add((relation_type, goid))

def reasoning(former, latter):
    if (former, latter) == ('is_a', 'is_a'):
        return 'is_a'
    elif (former, latter) in [('is_a', 'part_of'), ('part_of', 'is_a'), ('part_of', 'part_of')]:
        return 'part_of'
    elif (former, latter) in [('is_a', 'regulates'), ('regulates', 'is_a'), ('regulates', 'part_of')]:
        return 'regulates'
    else:
        return None

def get_inferred(inferred, former_rels, go_relations):
    for former_rel in former_rels:
        if go_relations.has_key(former_rel[1]):
            relations = set()
            for latter_rel in go_relations[former_rel[1]]:
                relation_type = reasoning(former_rel[0], latter_rel[0])
                if relation_type:
                    inferred.add(latter_rel[1])
                    relations.add((relation_type, latter_rel[1]))
            get_inferred(inferred, relations, go_relations)

def parse(handle):
    go_terms, go_inferreds = {}, {}

    go_term, go_relations = None, {}
    for line in handle:
        if line == '[Term]\n':
            go_term = GoTerm()
        elif go_term:
            if line.startswith('id: '):
                go_term.set_goid(line.split()[1])
            elif line.startswith('name: '):
                go_term.set_name(line[:-1].split(' ', 1)[1])
            elif line.startswith('is_a: '):
                go_term.add_relation('is_a', line.split()[1])
            elif line.startswith('relationship: part_of'):
                go_term.add_relation('part_of', line.split()[2])
            elif line.startswith('relationship: regulates') or \
                line.startswith('relationship: positively_regulates') or \
                line.startswith('relationship: negatively_regulates'):
                go_term.add_relation('regulates', line.split()[2])
            elif line.startswith('is_obsolete: '):
                go_term.obsolete = True
            elif line == '\n' and not go_term.obsolete:
                go_terms[go_term.goid] = go_term.name
                go_relations[go_term.goid] = go_term.relations
                go_term = None

    for (goid, relations) in go_relations.items():
        inferred = set([relation[1] for relation in relations])
        get_inferred(inferred, relations, go_relations)
        go_inferreds[goid] = inferred

    return go_terms, go_inferreds

if __name__ == '__main__':
    import sys
    from pprint import pprint

    go_terms, go_inferreds = parse(open(sys.argv[1]))

    print len(go_terms), len(go_inferreds)
    pprint(go_terms.items()[:5])
    pprint(go_inferreds.items()[:5])
