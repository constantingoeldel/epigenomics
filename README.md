# Extractor

This command line tool separates methylome files by the position of CG sites within gene body methylated genes. 

## Install from source: 

You need to have the rust compiler and the rust helper `cargo` installed (it's like npm to node).

On Linux, you can do a simple `curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh` to install both, for other OSs see https://www.rust-lang.org/learn/get-started

Clone the repository:

 `git clone git@github.com:constantingoeldel/extractor.git`

Compile and install the binary: 

`cargo install --path extractor`

To check it is successfully installed on your path, run `extractor -V`

## Run

Run `extractor --help` for always up-to-date argument syntax specification

![help options](help.png)

###  necessary arguments

Methylome directory: Path of directory containing the methlyome files from which to extract the CG-sites


`-m, --methylome <METHYLOME>`      


Annotation file:  Path of the annotation file containing information about beginning and end of gbM-genes

`-a, --annotation <ANNOTATION>   `

Output directory: Path of the directory where extracted segments shall be stored

`-o, --output-dir <OUTPUT_DIR>    `

### Optional arguments
Window size: Size of the window in percent of the gbM-gene length [default: 5]

 `-w, --window-size <WINDOW_SIZE>  `

Overwrite: Overwrite current content of the output directory?

 ` -f, --force `

Absolute size: Use absolute length in basis-pairs for window size instead of percentage of gene length?

 `--absolute                   `

Strandness: If supplied, will ignore the strand of the gene when determining where a give gene belongs

`-i, --ignore-strand`


## Examples: 

From `/mnt/extStorage/constantin/extractor` run 

```
extractor -m ../methylome/within_gbM_genes -a ../gbM_gene_anotation_extract_Arabidopsis.bed -w 5 -o ../windows 

extractor -m ../methylome/within_gbM_genes/ -g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed -o ../windows_relative_no_overlap  -w 5 -s 0 
```
# Development

For live-reloading changes use cargo watch 

```
cargo install cargo-watch
cargo watch -- cargo run -- -m ../methylome/within_gbM_genes/ -g ../methylome/gbM_gene_anotation_extract_Arabidopsis.bed -o ../windows_relative_no_overlap  -w 5 -s 0  
```

## Tests

For testing, run ```cargo test```

## Documentation

For documentation, run ```cargo doc --open```