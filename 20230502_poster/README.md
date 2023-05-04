# MITI 202305 Poster

- [Circular tree with background colored by phyla](data/generated/prune_tree/20230502_MITI_Genomes_Tree.circular_w_bgcolor.pruned.pdf)
- [Linear tree with background colored by phyla](data/generated/prune_tree/20230502_MITI_Genomes_Tree.rect_w_bgcolor.pruned.pdf)
- [Linear tree with no background](data/generated/prune_tree/20230502_MITI_Genomes_Tree.rect_no_color.pruned.pdf)

## Environment

1. Create a new virtual environment. Here is an example command using conda:

```bash
conda create -n trees python=3.11
```

2. Activate the environment:

```bash
conda activate trees
```

3. Install the required packages:

```bash
git clone git@github.com:FischbachLab/miti_small.git
cd miti_small/20230502_poster
pip install -U .
conda install -c conda-forge ncbi-datasets-cli
```

## How to generate the tree?

Once you have a list of species names that you'd like to visualize in a tree (example: [genomes.list](data/imported/genomes.list)). Upload them to the [NCBI Tax Identifier](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi) tool to get representative taxon ids. Then, create a file with one taxon id per line (see: [taxids.list](data/imported/taxids.list)). Finally, run the notebooks in the following order:

### Download the genomes

This notebook first identifies a representative genome for each taxon id provided as input, then downloads the genomes from NCBI using the `ncbi-datasets` command line tool.

Notebook: [download_genomes.ipynb](place_genome_on_tree/download_genomes.ipynb)

Data generated: `data/generated/download_genomes/`

```bash
data/generated/download_genomes/
├── genome_accessions.list
├── genomes.zip # (uploaded here: s3://genomics-workflow-core/Results/GTDB-tk-trees/MITI/20230502/genomes/)
└── genomes_metadata.csv
```

### Place genomes on the GTDB tree

We use gtdb-tk to place these genomes on a tree. Since this step requires a large amount of memory, I used a r5.12xlarge instance.

First, doenload the genomes from s3:

```bash
aws s3 sync s3://genomics-workflow-core/Results/GTDB-tk-trees/MITI/20230502/genomes genomes
```

Since we could not find a representative for the species *Harryflintia acetispora*, we downloaded the genome that we assembled from s3:

```bash
aws s3 cp s3://genomics-workflow-core/Results/HybridAssembly/MITI-MCB/SH0001499-00107/UNICYCLER/assembly.fasta genomes/Harryflintia_acetispora.MAF.fna
```

In hindsight, we could have done the same for *Blautia sp. (95% luti)* and *Blautia sp. (97% luti)*. But we can always rerun this if thats needed. Then, we run the [`gtdb_classify_wf.sh`](place_genome_on_tree/gtdb_scripts/gtdb_classify_wf.sh) script to place the genomes on the tree:

```bash
bash -x gtdb_classify_wf.sh $PWD/genomes
```

Data generated: `data/generated/process_genomes/`

```bash
data/generated/process_genomes
├── gtdb.20230502_MITI_Genomes_Tree.bac120.classify.tree
└── gtdb.20230502_MITI_Genomes_Tree.bac120.summary.tsv
```

### Prune the tree

This notebook uses the `treeViz.py` python script to prune the tree and generate a figure as a pdf file. It also procduces additional files (`.tree` - newick format tree) that can be used with a tree viewer such as [iTOL](https://itol.embl.de/) to make customizations. When using iTOL, you can use the [`display_names.csv`](data/imported/display_names.csv) file to map the taxon ids to the species names.

Notebook: [prune_tree.ipynb](place_genome_on_tree/prune_tree.ipynb)

Dependency: [treeViz.py](place_genome_on_tree/prune_tree.ipynb) python script

Data generated: `data/generated/prune_tree/`

```bash
data/generated/prune_tree/
├── 20230502_MITI_Genomes_Tree.circular_w_bgcolor.pruned.pdf
├── 20230502_MITI_Genomes_Tree.circular_w_bgcolor.pruned.tree
├── 20230502_MITI_Genomes_Tree.rect_no_color.pruned.pdf
├── 20230502_MITI_Genomes_Tree.rect_no_color.pruned.tree
├── 20230502_MITI_Genomes_Tree.rect_w_bgcolor.pruned.pdf
├── 20230502_MITI_Genomes_Tree.rect_w_bgcolor.pruned.tree
└── 20230502_MITI_Genomes_Tree.taxa_levels.csv
```

## Metadata

Google Sheet: `230502_MITI001_target_manufacturing.xlsx`