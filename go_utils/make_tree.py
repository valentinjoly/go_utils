#!/usr/bin/env python3

import argparse
import csv
import sys

from vjoly import (check_file_arg, open_file, print_stderr,
                   updict_append_to_list, STDIO_PATH)


PARENT_PREFIXES = [
    'is_a: GO',
    'intersection_of: GO',
    'intersection_of: part_of GO',
    'intersection_of: regulates GO',
    'intersection_of: positively_regulates GO',
    'intersection_of: negatively_regulates GO',
    'intersection_of: occurs_in GO',
    'relationship: part_of GO',
    'relationship: regulates GO',
    'relationship: positively_regulates GO',
    'relationship: negatively_regulates GO',
    'relationship: occurs_in GO',
]
CHILD_PREFIXES = [
    'intersection_of: has_part GO',
    'relationship: has_part GO',
]
INPUT_PATH = STDIO_PATH
OUTPUT_PATH = STDIO_PATH

TOP_GO_IDS = ('GO:0008150', 'GO:0003674', 'GO:0005575')


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)

    input_group = parser.add_argument_group(title='Input options')
    input_group.add_argument(
        'input_path', metavar='PATH', default=INPUT_PATH,
        help='Input OBO file. [{}]'.format(INPUT_PATH))

    output_group = parser.add_argument_group(title='Output options')
    output_group.add_argument(
        '-o', '--output_path', metavar='PATH', default=OUTPUT_PATH,
        help='Output GO tree file. [{}]'.format(OUTPUT_PATH))

    other_group = parser.add_argument_group(title='Other options')
    other_group.add_argument(
        '-h', '--help', action='help',
        help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    errors += check_file_arg(args.input_path, stdio_allowed=True)
    errors += check_file_arg(args.output_path, mode='w', stdio_allowed=True,
                             prefix='-o/--output_path')
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def import_obo(obo_path):
    terms, types, alt, alt_ids, obsolete = {}, {}, {}, {}, {}
    parents, children = {}, {}

    with open_file(obo_path, 'r') as obo_file:

        for line in obo_file:
            line = line.strip()

            if line == '[Typedef]':
                break

            elif line.startswith('id: '):
                go_id = line[4:]
                obsolete[go_id] = False
                alt[go_id] = False
                alt_ids[go_id] = []

            elif line.startswith('name: '):
                terms[go_id] = line[6:]
            elif line.startswith('namespace: '):
                types[go_id] = line[11:]
            elif line == 'is_obsolete: true':
                obsolete[go_id] = True

            elif line.startswith('alt_id: '):
                alt_id = line[8:]
                if '!' in alt_id:
                    alt_id = alt_id[:alt_id.index('!')]
                alt_id = alt_id.strip()
                alt[alt_id] = True
                updict_append_to_list(alt_ids, go_id, alt_id)

            else:
                for prefix in PARENT_PREFIXES:
                    if line.startswith(prefix):
                        parent = line[len(prefix) - 2:]
                        if '!' in parent:
                            parent = parent[:parent.index('!')]
                        parent = parent.strip()
                        updict_append_to_list(parents, go_id, parent)
                        updict_append_to_list(children, parent, go_id)

                for prefix in CHILD_PREFIXES:
                    if line.startswith(prefix):
                        child = line[len(prefix) - 2:]
                        if '!' in child:
                            child = child[:child.index('!')]
                        child = child.strip()
                        updict_append_to_list(parents, child, go_id)
                        updict_append_to_list(children, go_id, child)

    return terms, types, alt, alt_ids, obsolete, parents, children


def find_ancestors(parents):
    parents = {go_id: list(set(parents[go_id])) for go_id in parents}

    ancestors = {}
    for go_id in parents:
        prev_ancestors = []
        curr_ancestors = [go_id]

        while len(curr_ancestors) > len(prev_ancestors):
            prev_ancestors = curr_ancestors.copy()
            for ancestor in prev_ancestors:
                try:
                    curr_ancestors.extend(parents[ancestor])
                except KeyError:
                    continue
            curr_ancestors = list(set(curr_ancestors))

        curr_ancestors = list(set(curr_ancestors))
        curr_ancestors.pop(curr_ancestors.index(go_id))
        ancestors[go_id] = curr_ancestors

    return ancestors, parents


def eliminate_type_incongruities(types, parents, children, ancestors):
    for go_id, go_type in types.items():
        p = parents.get(go_id, None)
        if p is not None:
            parents[go_id] = [i for i in p if types[i] == go_type]

        c = children.get(go_id, None)
        if c is not None:
            children[go_id] = [i for i in c if types[i] == go_type]

        a = ancestors.get(go_id, None)
        if a is not None:
            ancestors[go_id] = [i for i in a if types[i] == go_type]


def assign_level(go_id, children, levels):
    try:
        for child in children[go_id]:
            try:
                levels[child] = max(levels[child], levels[go_id] + 1)
            except KeyError:
                levels[child] = levels[go_id] + 1
            assign_level(child, children, levels)
    except KeyError:
        return


def compute_go_levels(children, top_go_ids):
    levels = {}
    for top_go_id in top_go_ids:
        levels[top_go_id] = 1
        assign_level(top_go_id, children, levels)
    return levels


def add_alt_ids(terms, types, levels, obsolete, alt_ids,
                children, parents, ancestors):
    main_ids = {}
    for go_id in alt_ids:
        for alt_id in alt_ids[go_id]:
            main_ids[alt_id] = go_id
            terms[alt_id] = terms[go_id]
            types[alt_id] = types[go_id]
            obsolete[alt_id] = obsolete[go_id]
            try:
                levels[alt_id] = levels[go_id]
            except KeyError:
                pass
            for d in (children, parents, ancestors):
                try:
                    d[alt_id] = d[go_id].copy()
                except KeyError:
                    pass
    return main_ids


def export_tree(tree_path, terms, types, levels, obsolete, alt, main_ids,
                children, parents, ancestors):

    with open_file(tree_path, 'w') as tree_file:
        tree_table = csv.writer(tree_file, dialect='excel-tab')

        header = ["#GO_ID", "term", "type", "level", "is_obsolete",
                  "is_alternative", "main_id", "children", "parents",
                  "ancestors"]
        tree_table.writerow(header)

        for go_id in sorted(terms):
            try:
                level = levels[go_id]
            except KeyError:
                level = 'NA'
            try:
                main_id = main_ids[go_id]
            except KeyError:
                main_id = ''
            row = [
                go_id,
                terms[go_id],
                types[go_id],
                level,
                str(obsolete[go_id]),
                str(alt[go_id]),
                main_id
            ]
            for d in [children, parents, ancestors]:
                try:
                    row.append(', '.join(sorted(d[go_id])))
                except KeyError:
                    row.append('')
            tree_table.writerow(row)


def main(args):
    terms, types, alt, alt_ids, obsolete, parents, children = import_obo(
        args.input_path)

    ancestors, parents = find_ancestors(parents)
    eliminate_type_incongruities(types, parents, children, ancestors)
    levels = compute_go_levels(children, TOP_GO_IDS)

    main_ids = add_alt_ids(
        terms, types, levels, obsolete, alt_ids, children, parents, ancestors)

    export_tree(
        args.output_path, terms, types, levels, obsolete, alt, main_ids,
        children, parents, ancestors)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Make GO tree file.', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
