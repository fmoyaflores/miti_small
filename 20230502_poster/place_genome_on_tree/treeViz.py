#!/usr/bin/env python3
import argparse
from ete3 import PhyloTree, faces, AttrFace, TreeStyle, NodeStyle
import pandas as pd
import numpy as np
import seaborn as sns
import logging


def _get_colors(num_colors, color_palette="pastel"):
    """Generate N distinct colors.

    Args:
        num_colors (int): number of colors
        color_palette (str, optional): color (seaborn) palette to use. Defaults to 'pastel'.

    Returns:
        list: list of 'num_color' color hex values
    """
    return sns.color_palette(color_palette, num_colors).as_hex()


def parse_summary_file(summary_file, out_prefix, color=False, taxa_rank=None):
    """Parse the gtdb summary file to extract list of genomes and prune the tree,
    genomes per taxa rank and generate colors per taxa rank for the tree.

    Args:
        summary_file (str): path to summary file from GTDBtk classify_wf
        out_prefix (str): prefix for the parsed summary dataframe.
        taxa_rank ([type], optional): [description]. Defaults to phylum.

    Returns:
        tuple: containing 3 elements
            list:   list of all genomes in this subset
            dict:   taxa rank --> background color hex value
            dict:   list of genomes per taxa rank. taxa rank --> [genome1, genome2, ... , genomeN]
    """

    taxa_ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]

    if taxa_rank is None:
        taxa_rank = "phylum"

    taxa_rank = taxa_rank.lower()
    assert (
        taxa_rank in taxa_ranks
    ), f"Taxonomy rank must be one of [{','.join(taxa_ranks)}]"

    # Parse the summary file
    summary_df = pd.read_table(
        summary_file, header=0, usecols=["user_genome", "classification"]
    )
    summary_df[taxa_ranks] = summary_df["classification"].apply(
        lambda x: pd.Series(str(x).split(";"))
    )
    summary_df.to_csv(f"{out_prefix}.taxa_levels.csv", index=False)

    # Get the column number for the taxa rank for use in the itertuple for loop later
    taxa_rank_col_num = (
        taxa_ranks.index(taxa_rank) + 3
    )  # add 2 columns that already exist in the df

    genomes = sorted(summary_df["user_genome"].unique())
    logging.info(f"Final tree will contain {len(genomes)} genomes.")

    if color:
        color_dict = assign_color_by_rank(summary_df, taxa_rank)
        common_ancestor = get_common_ancestors(summary_df, taxa_rank_col_num)
    else:
        color_dict = None
        common_ancestor = None
    return genomes, color_dict, common_ancestor


def assign_color_by_rank(summary_df, taxa_rank):
    t_rank = sorted(summary_df[taxa_rank].unique())
    num_taxa_ranks = len(t_rank)

    colors = _get_colors(num_taxa_ranks) if num_taxa_ranks > 1 else ["#FFFFFF"]
    logging.info(f"Found {num_taxa_ranks} levels at {taxa_rank} rank.")

    # rank --> color
    color_dict = dict(zip(t_rank, colors))
    logging.info(color_dict)

    return color_dict


def get_common_ancestors(summary_df, taxa_rank_col_num):
    # rank --> [genome1, genome2, ... , genomeN]
    common_ancestor = {}
    for row in summary_df.itertuples():
        taxa = row[taxa_rank_col_num]
        if taxa in common_ancestor:
            common_ancestor[taxa].append(row.user_genome)
        else:
            common_ancestor[taxa] = [row.user_genome]
    return common_ancestor


def generate_tree(
    tree_file,
    genomes,
    out_prefix,
    circular=False,
    color_dict=None,
    common_ancestor=None,
    formatted_name_dict=None,
):
    """Read the tree from GTDBtk
      - prune it to the list of genomes provided.
      - color the nodes based on the taxa rank chosen
      - write a tree file in newick format and a draw a linear tree in pdf.

    Args:
        tree_file (str): large tree file from GTDBtk classify_wf.
        genomes (list): list of genomes to prune the large tree file.
        out_prefix (str): prexif for pruned tree file and pdf.
        circular (bool, optional): draw a circular tree. Defaults to False.
        color_dict (dict): taxa rank --> background color hex value
        common_ancestor (dict): list of genomes per taxa rank. taxa rank --> [genome1, genome2, ... , genomeN]
        fomatted_name_dict (dict): a dictionary of genome names and the names that should be printed in the figure
    """
    output_tree = f"{out_prefix}.pruned.tree"
    output_tree_fig_pdf = f"{out_prefix}.pruned.pdf"
    output_tree_fig_svg = f"{out_prefix}.pruned.svg"

    # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-backgrounds
    def layout(node):
        if node.is_leaf() and formatted_name_dict:
            full_name_face = faces.TextFace(formatted_name_dict[node.name])
            faces.add_face_to_node(full_name_face, node, column=0, position="aligned")

    # https://github.com/etetoolkit/ete/issues/194
    tree = PhyloTree(tree_file, format=1, quoted_node_names=True)

    logging.info(f"Pruning tree to {len(genomes)} organisms ...")
    tree.prune(genomes, preserve_branch_length=True)

    logging.info("Saving pruned tree in newick format ...")
    tree.write(format=1, outfile=output_tree)

    logging.info("Drawing pruned tree ...")
    ts = TreeStyle()
    ts.layout_fn = layout
    ts.show_leaf_name = False
    # ts.show_branch_length = True
    # ts.show_branch_support = True
    # ts.show_scale = True
    ts.complete_branch_lines_when_necessary = True
    ts.extra_branch_line_type = 0  # solid line
    ts.extra_branch_line_color = "black"
    ts.draw_guiding_lines = True
    ts.guiding_lines_type = 0  # solid line
    ts.guiding_lines_color = "black"
    ts.force_topology = True

    # Remove node glyphs from the tree
    remove_glyph = NodeStyle()
    remove_glyph["size"] = 0
    for d in tree.iter_descendants():
        d.img_style = remove_glyph

    ## Updating branch appearence
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#node-backgrounds
    if color_dict and common_ancestor:
        logging.info(f"Adding some color to the tree of {len(genomes)} organisms ...")
        for phyla, genome_list in common_ancestor.items():
            nstyle = NodeStyle()
            nstyle["bgcolor"] = color_dict[phyla]
            # nstyle["partition_bgcolor"] = color_dict[phyla]
            nstyle["hz_line_color"] = color_dict[phyla]
            nstyle["vt_line_color"] = color_dict[phyla]
            nstyle["hz_line_width"] = 1
            nstyle["vt_line_width"] = 1
            nstyle["size"] = 0
            if len(genome_list) == 1:
                tree.search_nodes(name=genome_list[0])[0].set_style(nstyle)
            else:
                tree.get_common_ancestor(genome_list).set_style(nstyle)
    else:
        nstyle = NodeStyle()
        nstyle["hz_line_color"] = "black"
        nstyle["vt_line_color"] = "black"
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_width"] = 1
        nstyle["size"] = 0

        for n in tree.traverse():
            n.set_style(nstyle)

    # For circular tree
    if circular:
        ts.mode = "c"
        ts.root_opening_factor = 1

    # tree.show(tree_style=ts)  # for interactive tree manipulation
    logging.info("Saving pruned tree image ...")
    tree.render(output_tree_fig_pdf, w=200, units="mm", tree_style=ts)
    return


def usage():
    parser = argparse.ArgumentParser(
        description="Generate a linear tree from GTDBtk tree file."
    )
    parser.add_argument(
        "-t",
        "--tree_file",
        help="Full tree file from GTDBtk classify_wf. NOTE: use the '--full_tree' option in GTDBtk v2.1.1 onwards.",
    )
    parser.add_argument(
        "-s",
        "--summary_file",
        help="Summary file from GTDBtk classify_wf. ",
    )
    parser.add_argument(
        "-o",
        "--out_prefix",
        help="Prefix for the output tree file and image.",
    )
    parser.add_argument(
        "-n",
        "--formatted_names",
        help="A csv file with two columns: node_name and display_name. The node_name should be the same as the genome name in the tree file. The display_name is the name that will be printed in the figure. If this option is not provided, the genome names will be printed in the figure.",
    )
    parser.add_argument(
        "--circular",
        action="store_true",
        help="Draw a circular tree.",
    )
    parser.add_argument(
        "--color_by_rank",
        action="store_true",
        help="Color the nodes based on the taxa rank. If true use the --highlight_rank option to highlight a specific taxa rank. Default: False.",
    )
    parser.add_argument(
        "--highlight_rank",
        default="phylum",
        help="Taxa rank to highlight. The default is 'phylum'.",
    )
    return parser.parse_args()


def main():
    args = usage()
    tree_file = args.tree_file
    summary_file = args.summary_file
    out_prefix = args.out_prefix
    formatted_names = args.formatted_names
    circular = args.circular
    color_by_rank = args.color_by_rank
    highlight_rank = args.highlight_rank

    if formatted_names is not None:
        fomatted_name_df = pd.read_csv(
            formatted_names,
            header=0,
            usecols=["node_name", "display_name"],
        )
        fomatted_name_dict = dict(
            zip(fomatted_name_df.node_name, fomatted_name_df.display_name)
        )
        logging.info("Tree leaves will be renamed according to the provided formatting")

    genomes, color_dict, common_ancestor = parse_summary_file(
        summary_file, out_prefix, color=color_by_rank, taxa_rank=highlight_rank
    )

    generate_tree(
        tree_file,
        genomes,
        out_prefix,
        circular,
        color_dict,
        common_ancestor,
        fomatted_name_dict,
    )
    logging.info("All done! Huzzah!")


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t[%(levelname)s]:\t%(message)s",
    )

    main()
