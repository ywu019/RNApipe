#!/usr/bin/env python

import re
import xml.etree.ElementTree as ET

def parse(handle, kobasdb):
    pathways, gene_pathways = [], set()

    tree = ET.parse(handle)
    root = tree.getroot()

    ncount, bcount = 0, 0
    for pathway in root.iter('Pathway'):
        source = pathway.find('Source').text
        if source == 'NCI-Nature Curated':
            ncount += 1
            pid = ncount * 10 + 2
            padb = 'n'
            paid = pathway.find('ShortName').text
        elif source == 'BioCarta':
            bcount += 1
            pid = bcount * 10 + 3
            padb = 'b'
            paid = pathway.get('id')
        else:
            continue

        pname = pathway.find('LongName').text
        pname = re.sub('&#038;#xef;', 'i', pname)
        pname = re.sub('&#xdf;', 'beta', pname)

        pathways.append((pid, padb, paid, pname))

        for pcomponent in pathway.iter('PathwayComponent'):
            interaction_id = pcomponent.get('interaction_idref')
            interaction_xpath = "./Model/InteractionList/Interaction[@id='" + interaction_id + "']"
            interaction = root.find(interaction_xpath)

            for icomponent in interaction.iter('InteractionComponent'):
                molecule_id = icomponent.get('molecule_idref')
                molecule_xpath = "./Model/MoleculeList/Molecule[@id='" + molecule_id + "']"
                molecule = root.find(molecule_xpath)
                molecule_type = molecule.get('molecule_type')

                gproteins = []  ##general proteins
                if molecule_type == 'protein':
                    gproteins.append(molecule)
                elif molecule_type == 'complex':
                    for ccomponent in molecule.iter('ComplexComponent'):
                        molecule_id = ccomponent.get('molecule_idref')
                        molecule_xpath = "./Model/MoleculeList/Molecule[@id='" + molecule_id + "']"
                        molecule = root.find(molecule_xpath)
                        molecule_type = molecule.get('molecule_type')
                        if molecule_type == 'protein':
                            gproteins.append(molecule)

                rproteins = []  ##real proteins, may be not true
                for gprotein in gproteins:
                    if gprotein.find('FamilyMemberList') is not None:
                        for member in gprotein.iter('Member'):
                            molecule_id = member.get('member_molecule_idref')
                            molecule_xpath = "./Model/MoleculeList/Molecule[@id='" + molecule_id + "']"
                            molecule = root.find(molecule_xpath)
                            rproteins.append(molecule)
                    elif gprotein.find('Part') is not None:
                        part = gprotein.find('Part')
                        molecule_id = part.get('whole_molecule_idref')
                        molecule_xpath = "./Model/MoleculeList/Molecule[@id='" + molecule_id + "']"
                        molecule = root.find(molecule_xpath)
                        rproteins.append(molecule)
                    else:
                        rproteins.append(gprotein)

                for rprotein in rproteins:
                    for name in rprotein.iter('Name'):
                        name_type = name.get('long_name_type')
                        if name_type == 'EntrezGene':
                            query = name.get('value')
                            gids = kobasdb.gids_from_entrez_gene_id(query)
                            for gid in gids:
                                gene_pathways.add((gid[0], pid))
                            break
                        elif name_type == 'UniProt':
                            query = name.get('value').split('-')[0]
                            gids = kobasdb.gids_from_uniprotkb_ac(query)
                            for gid in gids:
                                gene_pathways.add((gid[0], pid))
                            break

    pids = set()
    for gene_pathway in gene_pathways:
        pids.add(gene_pathway[1])

    for pathway in tuple(pathways):
        if pathway[0] not in pids:
            pathways.remove(pathway)

    pathways = tuple(pathways)
    gene_pathways = tuple(gene_pathways)

    return pathways, gene_pathways

if __name__ == '__main__':
    import sys
    from pprint import pprint

    from kobas import config, dbutils

    pathways, gene_pathways = parse(open(sys.argv[1]), \
        dbutils.KOBASDB(config.getrc()['kobasdb'] + 'hsa.db'))

    print len(pathways), len(gene_pathways)

    g, p = set(), set()
    for gp in gene_pathways:
        g.add(gp[0])
        p.add(gp[1])

    print len(g), len(p)

    pprint(pathways[:5])
    pprint(gene_pathways[:5])
