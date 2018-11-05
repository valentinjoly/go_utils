#!/usr/bin/env python3

import argparse
import csv
import sys

from vjoly import (check_file_arg, open_file, print_stderr,
                   updict_add_to_set, updict_set_union, STDIO_PATH)

INPUT_PATH = STDIO_PATH
OUTPUT_PATH = STDIO_PATH
MIN_LEVEL = 1
MAX_LEVEL = None


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, main=main)

    input_group = parser.add_argument_group(title='Input options')
    input_group.add_argument(
        'annot_paths', metavar='PATH', nargs='+', default=[INPUT_PATH],
        help='Input ANNOT file(s). [{}]'.format(INPUT_PATH))
    input_group.add_argument(
        '-t', '--tree_path', metavar='PATH', required=True,
        help='GO tree file.')
    input_group.add_argument(
        '-s', '--seqids_path', metavar='PATH',
        help='Table of sequence IDs to focus on.')

    annot_group = parser.add_argument_group(title='GO level selection')
    annot_group.add_argument(
        '-l', '--min_level', type=int, default=MIN_LEVEL,
        help='Minimum GO level to export annotations [{}]'.format(MIN_LEVEL))
    annot_group.add_argument(
        '-L', '--max_level', type=int, default=MAX_LEVEL,
        help='Maximum GO level to export annotations [{!s}]'.format(MAX_LEVEL))

    output_group = parser.add_argument_group(title='Output options')
    output_group.add_argument(
        '-o', '--output_path', metavar='PATH', default=OUTPUT_PATH,
        help='Output TSV table. [{}]'.format(OUTPUT_PATH))
    output_group.add_argument(
        '-a', '--add_ancestors', action='store_true',
        help='Add columns with ancestors of terminal GO terms.')

    other_group = parser.add_argument_group(title='Other options')
    other_group.add_argument(
        '-h', '--help', action='help',
        help='Print this help page and exit.')


def check_arguments(args):
    errors = []
    for annot_path in args.annot_paths:
        errors += check_file_arg(annot_path, stdio_allowed=True)
    errors += check_file_arg(args.tree_path,
                             prefix='-t/--tree_path')
    errors += check_file_arg(args.seqids_path, none_allowed=True,
                             prefix='-s/--seqids_path')
    errors += check_file_arg(args.output_path, mode='w', stdio_allowed=True,
                             prefix='-o/--output_path')
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)


def import_tree(tree_path):
    terms, types, levels = {}, {}, {}
    obsolete, main_ids, ancestors = {}, {}, {}
    with open_file(tree_path) as tree_file:
        for row in csv.reader(tree_file, dialect='excel-tab'):
            go_id = row[0]
            terms[go_id] = row[1]
            types[go_id] = row[2]
            try:
                levels[go_id] = int(row[3])
            except ValueError:
                levels[go_id] = None
            obsolete[go_id] = True if row[4] == 'True' else False
            if row[5] == 'True':
                main_ids[go_id] = row[6]
            ancestors[go_id] = set(row[-1].split(', '))
    return terms, types, levels, obsolete, main_ids, ancestors


def import_seqids(seqids_path):
    if seqids_path is not None:
        seqids = []
        with open_file(seqids_path) as seqids_file:
            for line in seqids_file:
                line = line.strip()
                seqids.append(line)
        return seqids
    return None


def add_annotation(annots, seqid, go_id, anc_annots, ancestors):
    updict_add_to_set(annots, seqid, go_id)
    if anc_annots is not None:
        updict_add_to_set(anc_annots, seqid, go_id)
        updict_set_union(anc_annots, seqid, ancestors[go_id])


def import_annotations(annot_paths, types, levels, obsolete, main_ids,
                       add_ancestors, ancestors, min_level=MIN_LEVEL,
                       max_level=MAX_LEVEL):
    bp, mf, cc, ec = {}, {}, {}, {}
    bp_anc, mf_anc, cc_anc = None, None, None
    if add_ancestors:
        bp_anc, mf_anc, cc_anc = {}, {}, {}

    for annot_path in annot_paths:
        with open_file(annot_path) as annot_file:
            for row in csv.reader(annot_file, dialect='excel-tab'):
                seqid = row[0]
                go_id = row[1]

                if not (go_id.startswith('GO:') or go_id.startswith('EC:')):
                    print_stderr(
                        'ERROR: Unknown annotation type: {}.'.format(go_id))
                    continue

                if go_id.startswith('EC:'):
                    updict_add_to_set(ec, seqid, go_id)
                    continue

                try:
                    go_type = types[go_id]
                except KeyError:
                    print_stderr(
                        'ERROR: Annotation not found: {}.'.format(go_id))
                    continue

                if obsolete[go_id]:
                    continue
                if levels[go_id] < min_level:
                    continue
                if max_level is not None and levels[go_id] > max_level:
                    continue
                try:
                    go_id = main_ids[go_id]
                except KeyError:
                    pass

                if go_type == 'biological_process':
                    add_annotation(bp, seqid, go_id, bp_anc, ancestors)
                elif go_type == 'molecular_function':
                    add_annotation(mf, seqid, go_id, mf_anc, ancestors)
                elif go_type == 'cellular_component':
                    add_annotation(cc, seqid, go_id, cc_anc, ancestors)

    return bp, mf, cc, ec, bp_anc, mf_anc, cc_anc


def format_terms(seqid, d, terms):
    try:
        go_ids = d[seqid]
    except KeyError:
        return ''
    text = []
    for go_id in sorted(go_ids):
        text.append('{} ({})'.format(terms[go_id], go_id))
    return '; '.join(text)


def format_ec_codes(seqid, ec):
    try:
        ec_codes = ec[seqid]
    except KeyError:
        return ''
    return '; '.join(sorted(ec_codes))


def export_table(output_path, seqids, terms, bp, mf, cc, ec,
                 add_ancestors, bp_anc, mf_anc, cc_anc):
    if seqids is None:
        seqids = sorted(list(set(bp.keys()) | set(mf.keys()) | set(cc.keys())))

    with open_file(output_path, 'w') as output_file:
        output_table = csv.writer(output_file, dialect='excel-tab')
        header = [
            'SeqID',
            'Biological process',
            'Molecular function',
            'Cellular component']
        if add_ancestors:
            header += [
                'Biological process (with ancestors)',
                'Molecular function (with ancestors)',
                'Cellular component (with ancestors)']
        header += ['Enzyme codes']
        output_table.writerow(header)
        for seqid in seqids:
            row = [
                seqid,
                format_terms(seqid, bp, terms),
                format_terms(seqid, mf, terms),
                format_terms(seqid, cc, terms)]
            if add_ancestors:
                row += [
                    format_terms(seqid, bp_anc, terms),
                    format_terms(seqid, mf_anc, terms),
                    format_terms(seqid, cc_anc, terms)]
            row.append(format_ec_codes(seqid, ec))
            output_table.writerow(row)


def main(args):
    terms, types, levels, obsolete, main_ids, ancestors = import_tree(
        args.tree_path)
    seqids = import_seqids(args.seqids_path)
    bp, mf, cc, ec, bp_anc, mf_anc, cc_anc = import_annotations(
        args.annot_paths, types, levels, obsolete, main_ids,
        args.add_ancestors, ancestors, args.min_level, args.max_level)
    export_table(args.output_path, seqids, terms, bp, mf, cc, ec,
                 args.add_ancestors, bp_anc, mf_anc, cc_anc)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert an ANNOT file to a per-sequence TSV table.',
        add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.main(args)
