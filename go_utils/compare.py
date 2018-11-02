#!/usr/bin/env python3

import argparse
import csv
import os
import scipy.stats
import sys

from go_utils.convert import import_tree

from vjoly import (
    check_bin_arg, check_dir_arg, check_file_arg, check_num_arg, open_file,
    print_stderr, run_child_process, updict_add_to_set, updict_append_to_list,
    updict_set_union, STDIO_PATH)


INPUT_PATH = STDIO_PATH
OUTPUT_DIR_PATH = '.'

TABLE = 'table'
TABLE_WITH_SEQIDS = 'table_with_seqids'
TABLE_WITH_SEQIDS_EXPANDED = 'table_with_seqids_expanded'
TEX = 'tex'
PDF = 'pdf'
OUTPUT_TYPE = TABLE
OUTPUT_TYPES = [OUTPUT_TYPE]
MIN_LEVEL = 1
MAX_LEVEL = None
MIN_FOLD_CHANGE = 1.5
MAX_PVALUE = 0.05
MIN_SEQ_COUNT = 1
MIN_SEQ_PROP = 0.0
LUALATEX_PATH = 'lualatex'

GO_TYPES = ('BP', 'MF', 'CC')
GO_TYPES_LONG = ('biological process',
                 'molecular function',
                 'cellular component')
TOP_GO_IDS = ('GO:0008150', 'GO:0003674', 'GO:0005575')

UP_REG = 'up'
DOWN_REG = 'down'
NOT_REG = 'not reg.'

TEX_HEADER = ["""\\documentclass[letterpaper,11pt]{article}

% LANGUAGE
\\usepackage{polyglossia}
\\setdefaultlanguage{english}

% FONT
\\usepackage{fontspec}
\\defaultfontfeatures{Ligatures=TeX,Scale=MatchLowercase}

% MARGINS
\\usepackage{geometry}
\\geometry{
  top=2cm,
  bottom=2cm,
  outer=2cm,
  inner=2cm
}
\\pagenumbering{gobble}
\\usepackage{setspace}

%TABLES
\\usepackage{float}
\\usepackage{tabularx}
\\usepackage{multirow}
\\usepackage{multicol}
\\usepackage{booktabs}
\\usepackage{longtable}
\\usepackage{ragged2e}
\\newcolumntype{P}[1]{>{\\RaggedRight\\hspace{0pt}}p{#1}}
\\usepackage{siunitx}

\\begin{document}
\\setlength\\parindent{0pt}
\\setlength\\tabcolsep{3pt}
\\setstretch{0.75}
\\renewcommand{\\arraystretch}{1.25}

GO enrichment analysis for “""", """” type annotations in sample “""", """”.
\\begin{footnotesize}
\\begin{longtable}[l]{@{}lP{5cm}ccrrcrrcrrc@{}}
\\toprule
  \\multirow{2}{*}[-2pt]{\\textbf{GO ID}}
& \\multirow{2}{*}[-2pt]{\\textbf{Description}}
& \\multirow{2}{*}[-2pt]{\\textbf{Level}}
&& \\multicolumn{2}{c}{\\textbf{Reference}}
&& \\multicolumn{2}{c}{\\textbf{""", """}}
&& \\multirow{2}{*}[-2pt]{\\textbf{FC}}
& \\multirow{2}{*}[-2pt]{\\textbf{\\emph{p}-value}}
& \\multirow{2}{*}[-2pt]{\\textbf{Reg.}}
\\\\
\\cmidrule{5-6}\\cmidrule{8-9}
&&&& \\multicolumn{1}{c}{\\textbf{Count}}
& \\multicolumn{1}{c}{\\textbf{Perc.}}
&& \\multicolumn{1}{c}{\\textbf{Count}}
& \\multicolumn{1}{c}{\\textbf{Perc.}}
&&&& \\\\
"""]

TEX_FOOTER = """\\bottomrule
\\end{longtable}
\\end{footnotesize}
\\end{document}
"""


def add_arguments(parser):
    parser.set_defaults(prog=parser.prog, check=check_arguments, main=main)

    input_group = parser.add_argument_group(title='Input parameters')
    input_group.add_argument(
        'annot_paths', metavar='PATH', nargs='+', default=[INPUT_PATH],
        help='Input ANNOT file(s). [{}]'.format(INPUT_PATH))
    input_group.add_argument(
        '-t', '--tree_path', metavar='PATH', required=True,
        help='GO tree file.')
    input_group.add_argument(
        '-r', '--ref_path', metavar='PATH', required=True,
        help='Path to list of reference sequence IDs.')
    input_group.add_argument(
        '-s', '--sample_path',  metavar='PATH', required=True,
        action='append', dest='sample_paths',
        help='Path to list of test sequence IDs.'
             '[repeat for multiple samples]')
    input_group.add_argument(
        '-S', '--sample_name', metavar='STR',
        action='append', dest='sample_names',
        help='Name of the test sample. [repeat for multiple samples]')
    input_group.add_argument(
        '-d', '--descriptions_path', metavar='PATH',
        help='Path to a table of sequence descriptions.')

    output_group = parser.add_argument_group(
        title='Output parameters [default type: {}]'.format(OUTPUT_TYPE))
    output_group.add_argument(
        '-o', '--output_dir_path', metavar='PATH', default=OUTPUT_DIR_PATH,
        help='Path to output directory. [{}]'.format(OUTPUT_DIR_PATH))
    output_group.add_argument(
        '--table', action='append_const',
        dest='output_types', const=TABLE,
        help='Normal table output.')
    output_group.add_argument(
        '--table_with_seqids', action='append_const',
        dest='output_types', const=TABLE_WITH_SEQIDS,
        help='Table output with an extra column for sequence IDs.')
    output_group.add_argument(
        '--table_with_seqids_expanded', action='append_const',
        dest='output_types', const=TABLE_WITH_SEQIDS_EXPANDED,
        help='Table output with one line per sequence ID, with descriptions.')
    output_group.add_argument(
        '--tex', action='append_const',
        dest='output_types', const=TEX,
        help='TeX output')
    output_group.add_argument(
        '--pdf', action='append_const',
        dest='output_types', const=PDF,
        help='TeX output, followed by compilation with LuaLaTeX')

    annot_group = parser.add_argument_group(title='Annotation filtering')
    annot_group.add_argument(
        '-P', '--process', action='store_true',
        help='Export "biological process" annotations.')
    annot_group.add_argument(
        '-F', '--function', action='store_true',
        help='Export "molecular function" annotations.')
    annot_group.add_argument(
        '-C', '--component', action='store_true',
        help='Export "cellular component" annotations.')
    annot_group.add_argument(
        '-l', '--min_level', type=int, default=MIN_LEVEL,
        help='Minimum GO level to export annotations [{}]'.format(MIN_LEVEL))
    annot_group.add_argument(
        '-L', '--max_level', type=int, default=MAX_LEVEL,
        help='Maximum GO level to export annotations [{!s}]'.format(MAX_LEVEL))

    enrichment_group = parser.add_argument_group(title='Enrichment criteria')
    enrichment_group.add_argument(
        '-f', '--min_fold_change', type=float, default=MIN_FOLD_CHANGE,
        help='Minimum absolute fold-change [{:.1f}]'.format(MIN_FOLD_CHANGE))
    enrichment_group.add_argument(
        '-p', '--max_pvalue', type=float, default=MAX_PVALUE,
        help='Maximum Fisher’s Exact Test p-value [{:.2f}]'.format(MAX_PVALUE))
    enrichment_group.add_argument(
        '-n', '--min_seq_count', type=int, default=MIN_SEQ_COUNT,
        help='Minimum sequence count [{:d}]'.format(MIN_SEQ_COUNT))
    enrichment_group.add_argument(
        '-m', '--min_seq_prop', type=float, default=MIN_SEQ_PROP,
        help='Minimum sequence proportion (%%) [{:.0f}]'.format(MIN_SEQ_PROP))

    regulation_group = parser.add_argument_group(
        title='GO terms to be reported [default: -UD]')
    regulation_group.add_argument(
        '-U', '--up_reg', action='store_true',
        help='Export up-regulated GO terms.')
    regulation_group.add_argument(
        '-D', '--down_reg', action='store_true',
        help='Export down-regulated GO terms.')
    regulation_group.add_argument(
        '-N', '--not_reg', action='store_true',
        help='Export not regulated GO terms.')

    other_group = parser.add_argument_group(title='Other parameters')
    other_group.add_argument(
        '-h', '--help', action='help',
        help='Print this help page and exit.')
    other_group.add_argument(
        '--lualatex_path', metavar='PATH', default=LUALATEX_PATH,
        help='Path to the LuaLaTeX executable [{}].'.format(LUALATEX_PATH))


def check_arguments(args):
    errors = []
    for annot_path in args.annot_paths:
        errors += check_file_arg(annot_path, stdio_allowed=True)
    for sample_path in args.sample_paths:
        errors += check_file_arg(sample_path, prefix='-s/--sample_path')
    errors += check_file_arg(args.tree_path,
                             prefix='-t/--tree_path')
    errors += check_file_arg(args.ref_path,
                             prefix='-r/--ref_path')
    errors += check_file_arg(args.descriptions_path, none_allowed=True,
                             prefix='-d/--descriptions_path')
    errors += check_dir_arg(args.output_dir_path, mode='w', create=True,
                            prefix='-o/--output_dir_path')
    errors += check_bin_arg(args.lualatex_path,
                            prefix='--lualatex_path')
    errors += check_num_arg(args.min_level, number_type=int,
                            mini=1, prefix='-l/--min_level')
    errors += check_num_arg(args.max_level, number_type=int,
                            mini=1, none_allowed=True, prefix='-L/--max_level')
    errors += check_num_arg(args.min_fold_change, number_type=float,
                            mini=1.0, prefix='-f/--min_fold_change')
    errors += check_num_arg(args.max_pvalue, number_type=float,
                            mini=0.0, maxi=1.0, prefix='-p/--max_pvalue')
    errors += check_num_arg(args.min_seq_count, number_type=int,
                            mini=0, prefix='-n/--min_seq_count')
    errors += check_num_arg(args.min_seq_prop, number_type=float,
                            mini=0.0, maxi=100.0, prefix='-m/--min_seq_prop')
    if args.max_level is not None and args.min_level > args.max_level:
        error = '-l/-L: Max level must be greater than min level.'
        errors.append(error)
    if args.sample_names is not None:
        if len(args.sample_names) != len(args.sample_paths):
            error = ('-s/-S: The number of sample names does not match '
                     'the number of sample paths.')
            errors.append(error)
    if errors:
        for error in errors:
            print_stderr(error, prefix='ERROR')
        sys.exit(1)

    if args.output_types is None:
        args.output_types = OUTPUT_TYPES

    if not (args.process or args.function or args.component):
        args.process, args.function, args.component = True, True, True

    if not (args.up_reg or args.down_reg or args.not_reg):
        args.up_reg, args.down_reg, args.not_reg = True, True, False


def get_name_from_path(path):
    name = os.path.basename(path)
    if '.' in name:
        return name[:name.rfind('.')]
    return name


def import_seqids(seqids_path):
    seqids = set()
    with open_file(seqids_path) as seqids_file:
        seqids_table = csv.reader(seqids_file, dialect='excel-tab')
        for row in seqids_table:
            if not row:
                continue
            seqids |= {row[0]}
    return seqids


def import_samples(sample_paths, sample_names):
    if sample_names is None:
        sample_names = [get_name_from_path(p) for p in sample_paths]
    samples = [import_seqids(p) for p in sample_paths]
    return samples, sample_names


def subtract_ref(samples, ref):
    return [ref - sample for sample in samples]


def add_annotation(annots, seqid, go_id, ancestors):
    updict_add_to_set(annots, seqid, go_id)
    updict_set_union(annots, seqid, ancestors[go_id])


def import_annotations(annot_paths, types, levels, obsolete, main_ids,
                       ancestors, min_level=MIN_LEVEL, max_level=MAX_LEVEL):
    bp, mf, cc = {}, {}, {}

    for annot_path in annot_paths:
        with open_file(annot_path) as annot_file:
            for row in csv.reader(annot_file, dialect='excel-tab'):
                seqid = row[0]
                go_id = row[1]

                if not go_id.startswith('GO:'):
                    if not go_id.startswith('EC:'):
                        print_stderr(
                            'ERROR: Unknown annotation type: '
                            '{}.'.format(go_id))
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
                    add_annotation(bp, seqid, go_id, ancestors)
                elif go_type == 'molecular_function':
                    add_annotation(mf, seqid, go_id, ancestors)
                elif go_type == 'cellular_component':
                    add_annotation(cc, seqid, go_id, ancestors)

    return bp, mf, cc


def import_descriptions(descriptions_path):
    if descriptions_path is None:
        return None
    descriptions = {}
    with open_file(descriptions_path) as descriptions_file:
        for row in csv.reader(descriptions_file, dialect='excel-tab'):
            seqid, desc = row
            descriptions[seqid] = desc
    return descriptions


def count_seqids(sample, annot):
    counts = {}
    for seqid in sample:
        seq_annot = annot.get(seqid, None)
        if seq_annot is None:
            continue
        for go_id in annot[seqid]:
            try:
                counts[go_id] += 1
            except KeyError:
                counts[go_id] = 1
    return counts


def add_zeroes(sample_counts, ref_counts):
    for go_id in sample_counts:
        ref_count = ref_counts.get(go_id, None)
        if ref_count is None:
            ref_counts[go_id] = 0
    for go_id in ref_counts:
        sample_count = sample_counts.get(go_id, None)
        if sample_count is None:
            sample_counts[go_id] = 0


def make_seqid_lists(sample, annot):
    seqids = {}
    for seqid in sample:
        seq_annot = annot.get(seqid, None)
        if seq_annot is None:
            continue
        for go_id in annot[seqid]:
            updict_append_to_list(seqids, go_id, seqid)
    return seqids


def compute_perc(count, top_count):
    if not top_count:
        return 0
    return 100 * count / top_count


def compute_fc(sample_perc, ref_perc):
    if not ref_perc:
        return float('+inf')
    elif not sample_perc:
        return float('-inf')
    elif sample_perc < ref_perc:
        return -1 * ref_perc / sample_perc
    return sample_perc / ref_perc


def compute_pvalue(sample_counts, ref_counts, go_id, top_go_id, pvalue_lookup):
    a = ref_counts[top_go_id]
    b = sample_counts[top_go_id]
    c = ref_counts[go_id]
    d = sample_counts[go_id]
    try:
        pvalue = pvalue_lookup[(a, b, c, d)]
    except KeyError:
        odds, pvalue = scipy.stats.fisher_exact([[a, b], [c, d]])
        pvalue_lookup[(a, b, c, d)] = pvalue
    return pvalue


def compute_reg(fc, pvalue, min_fc=MIN_FOLD_CHANGE, max_pvalue=MAX_PVALUE):
    if pvalue > max_pvalue:
        return NOT_REG
    if fc >= min_fc:
        return UP_REG
    if fc <= -1 * min_fc:
        return DOWN_REG
    return NOT_REG


def format_perc(perc):
    if not perc:
        return '0%'
    elif perc < 10:
        return '{:.2f}%'.format(perc)
    elif perc < 100:
        return '{:.1f}%'.format(perc)
    return '100%'


def format_fc(fc):
    return '{:+.3f}'.format(fc)


def format_pvalue(pvalue):
    # if pvalue < 0.0001:
    #     return '<0.0001'
    return '{:.4f}'.format(pvalue)


def process_sample(sample, ref, annot, top_go_id, terms, levels, descriptions,
                   min_seq_count, min_seq_prop, min_fc, max_pvalue,
                   export_not_reg, export_up_reg, export_down_reg,
                   pvalue_lookup):
    results = {}

    sample_counts = count_seqids(sample, annot)
    ref_counts = count_seqids(ref, annot)
    add_zeroes(sample_counts, ref_counts)
    seqids = make_seqid_lists(sample, annot)

    top_sample_count = sample_counts[top_go_id]
    top_ref_count = ref_counts[top_go_id]

    for go_id in sample_counts:
        term = terms[go_id]
        level = levels[go_id]
        sample_count = sample_counts[go_id]
        if sample_count < min_seq_count:
            continue
        ref_count = ref_counts[go_id]
        sample_perc = compute_perc(sample_count, top_sample_count)
        if sample_perc < min_seq_prop:
            continue
        ref_perc = compute_perc(ref_count, top_ref_count)
        if go_id == top_go_id:
            fc, pvalue, reg, key = '', '', '', ''
        else:
            fc = compute_fc(sample_perc, ref_perc)
            pvalue = compute_pvalue(
                sample_counts, ref_counts, go_id, top_go_id, pvalue_lookup)
            reg = compute_reg(fc, pvalue, min_fc, max_pvalue)
            if reg == NOT_REG and not export_not_reg:
                continue
            if reg == UP_REG and not export_up_reg:
                continue
            if reg == DOWN_REG and not export_down_reg:
                continue
            key = (-1*fc, go_id)
            fc = format_fc(fc)
            pvalue = format_pvalue(pvalue)
        sample_perc = format_perc(sample_perc)
        ref_perc = format_perc(ref_perc)
        result = [go_id, term, level, ref_count, ref_perc,
                  sample_count, sample_perc, fc, pvalue, reg]
        try:
            results[level][key] = result
        except KeyError:
            results[level] = {key: result}

    return results, seqids


def export_table(results, output_path):
    with open_file(output_path, 'w') as output_file:
        table = csv.writer(output_file, dialect='excel-tab')
        header = ['GO ID', 'Term', 'Level', 'Ref. count', 'Ref. perc.',
                  'Sample count', 'Sample perc.', 'FC', 'p-value', 'Reg.']
        table.writerow(header)
        for level in sorted(results):
            for key in sorted(results[level]):
                table.writerow(results[level][key])


def export_table_with_seqids(results, seqids, output_path):
    with open_file(output_path, 'w') as output_file:
        table = csv.writer(output_file, dialect='excel-tab')
        header = ['GO ID', 'Term', 'Level', 'Ref. count', 'Ref. perc.',
                  'Sample count', 'Sample perc.', 'FC', 'p-value', 'Reg.',
                  'Sequence IDs']
        table.writerow(header)
        for level in sorted(results):
            for key in sorted(results[level]):
                row = list(results[level][key])
                go_id = row[0]
                if level > 1:
                    row += [', '.join(sorted(seqids[go_id]))]
                table.writerow(row)


def export_table_with_seqids_expanded(results, seqids, descriptions,
                                      output_path):
    with open_file(output_path, 'w') as output_file:
        table = csv.writer(output_file, dialect='excel-tab')
        header = ['GO ID', 'Term', 'Level', 'Ref. count', 'Ref. perc.',
                  'Sample count', 'Sample perc.', 'FC', 'p-value', 'Reg.',
                  'Sequence ID']
        if descriptions is not None:
            header += ['Description']
        table.writerow(header)
        for level in sorted(results):
            for key in sorted(results[level]):
                row = results[level][key]
                go_id = row[0]

                lcl_seqids = seqids.get(go_id, None)
                if lcl_seqids is None or level < 2:
                    table.writerow(row)
                    continue
                lcl_seqids.sort()
                for i in range(len(lcl_seqids)):
                    if not i:
                        final_row = list(row)
                    else:
                        final_row = [''] * len(row)
                    final_row.append(lcl_seqids[i])
                    if descriptions is not None:
                        desc = descriptions.get(lcl_seqids[i], None)
                        if desc is not None:
                            final_row.append(desc)
                    table.writerow(final_row)


def format_tex_row(row):
    row[0] = '\\textbf{' + row[0] + '}'
    for i in (3, 6, 9):
        row.insert(i, '')
    if row[10] == '+inf':
        row[10] = '$+\\infty$'
    if row[10] == '-inf':
        row[10] = '$-\\infty$'
    else:
        row[10] = '$' + row[10] + '$'
    row[11] = '$' + row[11] + '$'
    row = ' & '.join([str(i) for i in row]) + ' \\\\'
    row = row.replace(' &  & ', ' && ')
    row = row.replace('%', '\\%')
    row = row.replace('_', '\\_')
    return row


def export_tex(results, output_path, name, go_type_long,
               header=TEX_HEADER, footer=TEX_FOOTER):
    name = name.replace('_', '\\_')
    go_type_long = go_type_long.replace('_', '\\_')
    with open_file(output_path, 'w') as output_file:
        output_file.write(''.join(
            [header[0], go_type_long, header[1], name,
             header[2], name, header[3]]))
        for level in sorted(results):
            output_file.write('\\midrule\n')
            for key in sorted(results[level]):
                row = format_tex_row(results[level][key])
                output_file.write(row + '\n')
        output_file.write(footer)


def compile_tex(output_path, output_dir_path=OUTPUT_DIR_PATH,
                lualatex_path=LUALATEX_PATH):
    lualatex_cmd = [lualatex_path,
                    '--output-directory=' + output_dir_path,
                    output_path]
    run_child_process(lualatex_cmd, remove_log=True)


def export_results(results, name, go_type, go_type_long, seqids, descriptions,
                   output_dir_path=OUTPUT_DIR_PATH, output_types=OUTPUT_TYPES,
                   lualatex_path=LUALATEX_PATH):
    output_basename = '_'.join(['GO', 'enrichment', name, go_type])
    output_path_base = os.path.join(output_dir_path, output_basename)
    if TABLE in output_types:
        output_path = output_path_base + '.txt'
        export_table(results, output_path)
    if TABLE_WITH_SEQIDS in output_types:
        output_path = output_path_base + '+seqids.txt'
        export_table_with_seqids(results, seqids, output_path)
    if TABLE_WITH_SEQIDS_EXPANDED in output_types:
        output_path = output_path_base + '+seqids_exp.txt'
        export_table_with_seqids_expanded(
            results, seqids, descriptions, output_path)
    if TEX in output_types or PDF in output_types:
        output_path = output_path_base + '.tex'
        export_tex(results, output_path, name, go_type_long)
        if PDF in output_types:
            compile_tex(output_path, output_dir_path, lualatex_path)


def main(args):
    tree = import_tree(args.tree_path)
    terms, types, levels, obsolete, main_ids, ancestors = tree

    annotations = import_annotations(
        args.annot_paths, types, levels, obsolete, main_ids,
        ancestors, args.min_level, args.max_level)

    samples, names = import_samples(args.sample_paths, args.sample_names)
    descriptions = import_descriptions(args.descriptions_path)

    ref = import_seqids(args.ref_path)
    refs = subtract_ref(samples, ref)

    pvalue_lookup = {}
    analyze_types = (args.process, args.function, args.component)

    for analyze_type, annot, go_type, go_type_long, top_go_id in zip(
            analyze_types, annotations, GO_TYPES, GO_TYPES_LONG, TOP_GO_IDS):
        if not analyze_type:
            continue
        for sample, ref, name in zip(samples, refs, names):
            results, seqids = process_sample(
                sample, ref, annot, top_go_id, terms, levels,
                descriptions, args.min_seq_count, args.min_seq_prop,
                args.min_fold_change, args.max_pvalue, args.not_reg,
                args.up_reg, args.down_reg, pvalue_lookup)

            export_results(
                results, name, go_type, go_type_long, seqids, descriptions,
                args.output_dir_path, args.output_types, args.lualatex_path)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Perform GO term enrichment analyses.',
        add_help=False)
    add_arguments(parser)
    args = parser.parse_args()
    args.check(args)
    args.main(args)
