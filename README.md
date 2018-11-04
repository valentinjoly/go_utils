# GO utils

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
  - [Build the internal GO tree with `make_tree`](build-the-internal-go-tree-with-make_tree)
  - [Convert annotation files with `convert`](#convert-annotation-files-with-convert)
  - [Perform enrichment analyses with `compare`](#perform-enrichment-analyses-with-compare)
- [License](#license)

## Introduction
**GO utils** is a small set of GO (*Gene Ontology*) analysis utilities allowing to manipulate functional annotation data. It comprises three different programs:
- **`make_tree`** builds an internal GO annotation tree used as an input for `convert` and `compare`;
- **`convert`** converts an `.annot` file produced by your annotation program (e.g. Blast2GO) into a tab-separated table;
- **`compare`** performs GO-term enrichment analyses between a sample and a reference.

## Installation
GO utils is coded in **Python 3**.  Check out *The Hitchhiker’s Guide to Python* for installation prpocedures on [Linux](https://docs.python-guide.org/starting/install3/linux/#install3-linux), [Windows](https://docs.python-guide.org/starting/install3/win/#install3-windows), and [Mac OS X](https://docs.python-guide.org/starting/install3/osx/#install3-osx).

GO utils also requires **SciPy**. Detailed installation instructions are available on the [SciPy website](https://www.scipy.org/install.html).

GO utils relies on some functions from my own all-purpose module, **`vjoly.py`**, available from my GitHub [`vjoly`](https://github.com/valentinjoly/vjoly) repository. Clone or download this repo and add `vjoly.py` to your `PYTHONPATH` environment variable. For instance, on Linux:

    $ cd <installation_path>
    $ git clone git@github.com:valentinjoly/vjoly.git
    $ echo 'export PYTHONPATH=<installation_path>/vjoly:$PYTHONPATH' >> ~/.bashrc
    $ source ~/.bashrc

Finally, clone or download the **`go_utils`** repository and add the `go_utils/bin` subdirectory to your `PATH` variable. For instance, on Linux:

    $ cd <installation_path>
    $ git clone git@github.com:valentinjoly/go_utils.git
    $ echo 'export PATH=<installation_path>/go_utils/bin:$PATH' >> ~/.bashrc
    $ source ~/.bashrc

## Usage

### Build the internal GO tree with `make_tree`

**`make_tree`** builds an internal GO annotation tree used as an input for `convert` and `compare`. The input is the core GO ontology file in OBO format provided by the [Gene Ontology Consortium](http://geneontology.org/page/download-ontology). The permanent link to download the **`go.obo`** file is: [http://purl.obolibrary.org/obo/go.obo](http://purl.obolibrary.org/obo/go.obo).

To build the internal `go_utils` tree, use one of the following commands:

    $ go_utils make_tree go.obo > go_tree.txt
    $ go_utils make_tree -o go_tree.txt go.obo
    $ cat go.obo | go_utils make_tree - > go_tree.txt

### Convert annotation files with `convert`

#### Basic usage

**`convert`** converts an `.annot` file produced by your annotation program (e.g. Blast2GO) into a tab-separated table.

The minimum command for `convert` is:

    $ go_utils convert -t go_tree.txt [options] annotations.annot > annotations_table.txt

Several `.annot` files can be specified at once in the command line:

    $ go_utils convert -t go_tree.txt [options] annot1.annot annot2.annot annot3.annot > annotations_table.txt

In all cases, a unique tab-separated table will be produced (file `annotations_table.txt` in the example above). The table looks like:

| **SeqID** | **Biological process** | **Molecular function** | **Cellular component** | **Enzyme codes** |
| --- | --- | --- | --- | --- |
| **seq1** | go_term_1 (GO:0000001); go_term_2 (GO:0000002) | go_term_3 (GO:0000003); go_term_4 (GO:0000004); go_term_5 (GO:0000005) | | EC:0.0.0.0 |
| **seq2** | go_term_1 (GO:0000001); go_term_6 (GO:0000006); go_term_7 (GO:0000007) | go_term_3 (GO:0000003) | go_term_8 (GO:0000008); go_term_9 (GO:0000009) | |
| **seq3** | go_term_2 (GO:0000002); go_term_10 (GO:0000010) | go_term_11 (GO:0000011) | go_term_12 (GO:0000012) | EC:0.0.0.1; EC:0.0.0.2 |

The following adjustments are  made:
- GO terms tagged as *obsolete* in the initial `go.obo` file are ignored;
- GO terms tagged as *alternative* in the initial `go.obo` file are replaced with the appropriate representative GO term;
- Ancestors of a GO term that do not belong to the same GO type (*biological process*, *molecular function*, *cellular component*) are ignored if option `-a` is set.

#### Optional parameters

**`-o/--output_path`** allows to specify the path to the file where `convert` will write the output table. By default, `convert` writes on the standard output.

**`-s/--seqids_path`** allows to specify the path to a text file containing a list of sequence IDs of interest (one ID per line). By default, `convert` export annotations for all sequence IDs present in the input `.annot` file(s). Option `-s` allows to narrow down on a smaller set of sequences. Also, if sequences found in this list are not found in the `.annot` files(s), `convert` will add empty rows for them in the output table.

**`-a/--add_ancestors`** allows to add three additional columns in the output table, where `convert` will write not only the terminal GO terms found in the `.annot` file(s), but also their ancestors in the GO tree.

**`-l/--min_level`** and **`-L/--max_level`** allow to specify the minimum and maximum hierarchical GO levels for exported annotations. By default, `convert` exports annotations from all levels (`min_level = 1`, `max_level = None`)


### Perform enrichment analyses with `compare`

**`compare`** performs GO-term enrichment analyses between one or more samples and a reference.

#### Algorithm details

1. Import the `*.annot` file(s) to determine which terminal GO terms are associated to each sequence.
2. Use the internal GO tree to determine which ancestor GO terms are associated to each sequence.
3. Count the number of sequences from the sample of interest and from the reference that are associated to each GO term.
4. For each GO-term *x*, compute the following values:
  - **_S<sub>x</sub>_** : the number of sequences associated to GO term *x* in the sample
  - **_R<sub>x</sub>_** : the number of sequences associated to GO term *x* in the reference
  - **_S_<sub>top</sub>** : the number of sequences associated to the corresponding top GO term* in the sample
  - **_R_<sub>top</sub>** : the number of sequences associated to the corresponding top GO term* in the reference
  - **_P<sub>x</sub>_ = _S<sub>x</sub>_ / _S_<sub>top</sub>** : the proportion of annotated sequences associated to GO term *x* in the sample
  - **_Q<sub>x</sub>_ = _R<sub>x</sub>_ / _R_<sub>top</sub>** : the proportion of annotated sequences associated to GO term *x* in the reference
5. Compute the enrichment values:
  - **_D<sub>x</sub>_ = _P<sub>x</sub>_ / _Q<sub>x</sub>_** : the enrichment ratio between the sample and the reference
  - **_FC<sub>x</sub>_** : the corresponding fold-change, given by *FC<sub>x</sub>* = *D<sub>x</sub>* if *D<sub>x</sub>* ≥ 1 and *FC<sub>x</sub>* = −1/*D<sub>x</sub>* if *D<sub>x</sub>* < 1
6. Assess the statistical significance of the enrichment by computing the **_p_-value** of a Fisher’s Exact Test on the [(_S<sub>x</sub>_, _R<sub>x</sub>_), (_S_<sub>top</sub>, _R_<sub>top</sub>)] contigency matrix.
7. Determine whether the GO-term is enriched or not. Up-regulated GO-terms are defined by **_FC_ ≥ 1.5** and **_p_ ≤ 0.05**; down-regulated GO-terms are defined by  **_FC_ ≤ −1.5** and **_p_ ≤ 0.05**. Default selection criteria can be modified with options `-f` and `-p`.
8. Write results for each sample and each GO type in a distinct table. Different output formats are available (see below).

*Top GO terms are: *biological_process* (GO:0008150), *molecular_function* (GO:0003674), and *cellular component* (GO:0005575).

The default output is a tab-delimited table that looks like:

| GO ID      | Term               | Level | Ref. count | Ref. perc. | Sample count | Sample perc. | FC  | p-value | Reg. |
| ---------- | ------------------ | ----- | ---------- | ---------- | ------------ | ------------ | --- | ------- | -----|
| GO:0008150 | biological_process | 1 | *D* | 100% | *C* | 100% | | | |
| GO:0000001 | go_description_1   | 2 | *B*<sub>1</sub> |*Q*<sub>1</sub> | *A*<sub>1</sub> | *P*<sub>1</sub> | *FC*<sub>1</sub> | *p*<sub>1</sub> | UP   |
| GO:0000002 | go_description_2   | 2 | *B*<sub>2</sub> |*Q*<sub>2</sub> | *A*<sub>2</sub> | *P*<sub>2</sub> | *FC*<sub>2</sub> | *p*<sub>2</sub> | DOWN |
| GO:0000003 | go_description_3   | 3 | *B*<sub>3</sub> |*Q*<sub>3</sub> | *A*<sub>3</sub> | *P*<sub>3</sub> | *FC*<sub>3</sub> | *p*<sub>3</sub> | UP   |

As with `convert`, the following adjustments are made:
- GO terms tagged as *obsolete* in the initial `go.obo` file are ignored;
- GO terms tagged as *alternative* in the initial `go.obo` file are replaced with the appropriate representative GO term;
- Ancestors of a GO term that do not belong to the same GO type (*biological process*, *molecular function*, *cellular component*) are ignored.

#### Basic usage

The minimum command for `compare` is:

    $ go_utils compare -t go_tree.txt -s sample.txt -r reference.txt annotations.annot

**`-t/--tree_path`** : path to the GO tree built with `make_tree`.

**`-r/--ref_path`** : path to the list of sequence IDs constituting the reference for the enrichment analysis (one ID per line).

**`-s/--sample_path`** : path to the list of sequence IDs constituting the sample of interest for the enrichment analysis (one ID per line). Option `-s` can be used multiple times to specify different samples of interest. A distinct analysis will be performed for each of them against the same reference.

For each sample, the final set of reference sequences for the enrichment analysis is defined as the sequence set specified with option `-r` minus the sequence set specified with option `-s`.

Several `.annot` files can be specified at once in the command line:

    $ go_utils compare -t go_tree.txt -s sample.txt -r reference.txt annot1.annot annot2.annot annot3.annot

#### Optional parameters

**`-S/--sample_name`** can be used to give a name to the sample specified with option `-s`. By default, `compare` deduces the sample name from the file name given with option `-s`. If several samples were specified with option `-s`, option `-S` must be repeated the same number of times, in the same order.

**`-i/--info_path`** allows to specify additional information about each sequence in tab-separated format, with a mandatory header. This information will be displayed in the output table if option `--table_with_seqids_expanded` is set.

**`-o/--output_dir_path`** allows to specify the path to the directory where `compare` will write the output file(s). It can be an existing directory, or the path to a new directory which will be created. By default, `convert` writes in the current working directory.

##### Enrichment criteria

**`-f/--min_fold_change`** and **`-p/--max_pvalue`** can be used to modify the minimum absolute fold-change and maximum Fisher’s exact test *p*-value, respectively, defining enriched GO-terms. The default is `-f 1.5 -p 0.05`. The minimum fold-change must be a real number ≥ 1. The maximum *p*-value must be a real number between 0 and 1 (inclusive).

**`-n/--min_seq_count`** and **`-m/--min_seq_prop`** allow to specify the minimum number of sequences and the minimum proportion of sequences, respectively, associated to a GO term to consider it for enrichment analyses. If both options are set, a GO term must satisfy both criteria to be selected. The default is `-n 1 -m 0`.

##### Selection of exported GO-terms

**`-P/--process`**, **`-F/--function`**, and **`-C/--component`** can be used to specify on which type(s) of GO annotations `compare` should perform enrichment analyses. If more than one of these options are set, `convert` will process each GO type separately. Setting none of these options is equivalent to using `-PFC`.

**`-l/--min_level`** and **`-L/--max_level`** allow to specify the minimum and maximum hierarchical GO levels for enrichment analyses. By default, `compare` analyzes annotations from all levels (`min_level = 1`, `max_level = None`)

**`-U/--up_reg`**, **`-D/--down_reg`**, and **`-N/--not_reg`** allow to export up-, down-, or not regulated GO terms in the output. By default, `compare` exports up- and down-regulated GO-terms. (`-UD`)

##### Output formats

One or several of the options below can be used to specify output format(s). The default is `--table`.

**`--table`** is the default option, for a tab-separated output as shown above.

**`--table_with_seqids`** produces the same output as `--table` with an extra column on the right containing a comma-separated list of sequence IDs associated to each enriched GO term. By default, this column is populated for enriched GO-terms only. However, option **`--export_level_one_seqids`** can be set if one wants to get the list of sequence IDs associated to level 1 GO-terms (i.e., `biological_process`, `molecular_function`, and `cellular_component`).

**`--table_with_seqids_expanded`** produces the same output as `--table_with_seqids`, except that the list of sequence IDs associated to each GO term is distributed across several rows, each of which can be optionally associated to tabular sequence information provided with option `--info_path`. Here again, option **`--export_level_one_seqids`** can be set if one wants to include the list of sequence IDs associated to level 1 GO-terms (i.e., `biological_process`, `molecular_function`, and `cellular_component`).

**`--tex`** produces the same output as `--table`, but formatted in TeX. Compilation of this file requires LuaLaTeX and the following packages: `booktabs`, `float`, `fontspec`, `geometry`, `longtable`, `multicol`, `multirow`, `polyglossia`, `ragged2e`, `setspace`, `siunitx`, and `tabularx`.

**`--pdf`** produces the same TeX file as with option `--tex` and compiles it with LuaLaTeX to make a PDF document. The `lualatex` program and appropriate packages (see above) must be installed on your system. Unless specified with option **`--lualatex_path`**, path to the `lualatex` executable is assumed to be part of your `PATH` variable. 

## License
This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE](LICENSE) file for details.
