from Bio import Entrez
import pandas as pd
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Align.Applications import ClustalwCommandline
import re
import time


def alignment(dire, file, output_file):
    """
    Executes a multi-sequence alignment using ClustalW.
    :param dire: diretory to the ClustalW tool.
    :param file: list of fasta files to perform the alignment.
    :param output_file: str with the name to the output file.
    :return:
    """
    clustalw_cline = ClustalwCommandline(dire, infile=file)
    clustalw_cline()
    print(clustalw_cline)
    print(output_file + '.aln e ' + output_file + '.dnd gerados.')



# diretoria = r"/path/to/tool/clustalw2"
#
# dir_file = [r"/path/to/file/Fasta_baseplate.fasta", r"/path/to/file/Fasta_tail.fasta", r"/path/to/file/Fasta_tail_fiber.fasta", r"/path/to/file/Fasta_tail_sheath.fasta"]
# out_file = ["Baseplate_aln", "Tail_aln", "Tail_fiber_aln", "Tail_sheath_aln"]


# for i in range(len(dir_file)):
#     alignment(diretoria, dir_file[i], output_file=out_file[i])

def read_aligment(file):
    """
    Reads the alignemnt performed with ClustalW.
    :param file: list of alignment files to read.
    :return:
    """
    align = AlignIO.read(file, "clustal")
    print(align)

def phylogeny(file):
    """
    Prints the phylogenetic tree in ascci.
    :param file: list of files to draw the trees.
    :return:
    """
    tree = Phylo.read(file, "newick")
    Phylo.draw_ascii(tree)

#
# aln_files = ["Fasta_tail.aln", "Fasta_tail_fiber.aln", "Fasta_tail_sheath.aln"]
# dnd_files = ["Fasta_baseplate.dnd", "Fasta_tail.dnd", "Fasta_tail_fiber.dnd", "Fasta_tail_sheath.dnd"]

#
# start_time = time.time()

# for i in range(len(aln_files)):
#     read_aligment(aln_files[i])

# for j in range(len(dnd_files)):
#     phylogeny(dnd_files[j])
#
# if __name__ == "__main__":
#     print("--- %s seconds ---" % (time.time() - start_time))


