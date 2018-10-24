#!/usr/bin/env python3

import argparse
import csv
import sys

from vjoly import check_file_arg, open_file, print_stderr, STDIO_PATH


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


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, process=main)
    parser.add_argument(
        'input_path', metavar='FILE', default=INPUT_PATH,
        help='Input OBO file. [{}]'.format(INPUT_PATH))
    parser.add_argument(
        '-o', '--output_path', metavar='FILE', default=OUTPUT_PATH,
        help='Output GO tree file. [{}]'.format(OUTPUT_PATH))
    parser.add_argument(
        '-h', '--help', action='help',
        help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    errors += check_file_arg(args.input_path, stdio_allowed=True)
    errors += check_file_arg(args.output_path, mode='w', stdio_allowed=True)
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def import_obo(obo_path):

    terms, types, alt, alt_ids, obsolete, parents = {}, {}, {}, {}, {}, {}

    with open_file(obo_path, 'r') as obo_file:

        for line in obo_file:
            line = line.strip()

            if line == '[Typedef]':
                break

            elif line.startswith('id: '):
                go_id = line[4:]
                obsolete[go_id] = False
                alt[go_id] = False
                parents[go_id] = []
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
                alt_ids[go_id].append(alt_id)
                alt[alt_id] = True

            else:
                for prefix in PARENT_PREFIXES:
                    if line.startswith(prefix):
                        go_parent = line[len(prefix) - 2:]
                        if '!' in go_parent:
                            go_parent = go_parent[:go_parent.index('!')]
                        go_parent = go_parent.strip()
                        parents[go_id].append(go_parent)

                for prefix in CHILD_PREFIXES:
                    if line.startswith(prefix):
                        go_child = line[len(prefix) - 2:]
                        if '!' in go_child:
                            go_child = go_child[:go_child.index('!')]
                        go_child = go_child.strip()
                        if go_child not in parents:
                            parents[go_child] = []
                        parents[go_child].append(go_id)

    return terms, types, alt, alt_ids, obsolete, parents


def find_ancestors(parents):

    parents = {go_id: list(set(parents[go_id])) for go_id in parents}

    ancestors = {}
    for go_id in parents:
        prev_ancestors = []
        curr_ancestors = [go_id]

        while len(curr_ancestors) > len(prev_ancestors):
            prev_ancestors = curr_ancestors.copy()
            for ancestor in prev_ancestors:
                curr_ancestors.extend(parents[ancestor])
            curr_ancestors = list(set(curr_ancestors))

        curr_ancestors = list(set(curr_ancestors))
        curr_ancestors.pop(curr_ancestors.index(go_id))
        ancestors[go_id] = curr_ancestors

    return ancestors, parents


def add_alt_ids(terms, types, alt_ids, obsolete, parents,
                ancestors):

    for go_id in alt_ids:
        for alt_id in alt_ids[go_id]:
            terms[alt_id] = terms[go_id]
            types[alt_id] = types[go_id]
            obsolete[alt_id] = obsolete[go_id]
            parents[alt_id] = parents[go_id].copy()
            ancestors[alt_id] = ancestors[go_id].copy()

    return terms, types, obsolete, parents, ancestors


def export_tree(tree_path, terms, types, alt, obsolete, parents, ancestors):

    with open_file(tree_path, 'w') as tree_file:
        tree_table = csv.writer(tree_file, dialect='excel-tab')

        header = ["#GO_ID", "term", "type", "is_alternative", "is_obsolete",
                  "parents", "ancestors"]
        tree_table.writerow(header)

        for go_id in sorted(terms):
            row = [
                go_id,
                terms[go_id],
                types[go_id],
                str(alt[go_id]),
                str(obsolete[go_id]),
                ', '.join(sorted(parents[go_id])),
                ', '.join(sorted(ancestors[go_id]))
            ]
            tree_table.writerow(row)


def main(args):

    terms, types, alt, alt_ids, obsolete, parents = import_obo(args.input_path)
    ancestors, parents = find_ancestors(parents)
    terms, types, obsolete, parents, ancestors = add_alt_ids(
        terms, types, alt_ids, obsolete, parents, ancestors)
    export_tree(args.output_path, terms, types, alt, obsolete,
                parents, ancestors)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Make GO tree file', add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.process(args)
