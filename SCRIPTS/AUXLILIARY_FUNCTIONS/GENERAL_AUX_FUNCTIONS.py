from hashlib import new
import re
from nbformat import write
import pandas as pd
from Bio import SeqFeature
import os


def remove_id(list_ids, unwanted_ids):
    """From a given list of ids remove the the ids inside the list_ids that have the expressions prsented in the unwanted_list.

    Args:
        list_ids (list): list containing all the ids
        unwanted_ids (list): list of str containing the unwanted expressions in the ids to be removed.

    Returns:
    """
    remove_id = []
    for i in list_ids:
        if unwanted_ids in i:
            remove_id.append(i)
    for j in remove_id:
        if j in list_ids:
            list_ids.remove(j)
    return list_ids



def sort_by_id(list_ids, ids):
    new = []
    for i in list_ids:
        if ids in i:
            new.insert(0, i)
        if ids not in i:
            new.append(i)
    return new


def get_feature(seqio_obje, feature):
    """
    Allows to retrieve features from a genbank file.
    :param seqio_obje: object with the information of the genbank file
    :param feature: type of feature to extract
    :return: the wanted feature
    """
    for feat in seqio_obje.features:
        if feat.type == feature:
            feat_type = feat
            return feat_type


def remove_string_prefix(string, prefix):
    if string.startswith(prefix):
        return string[len(prefix):]
    else:
        return None


def remove_string_suffix(string, suffix):
    if string.endswith(suffix):
        print(len(suffix))
        return string[:len(suffix)]
    else:
        return None


def id_index_in_dataset(in_csv, ids, column):
    df = pd.read_csv(in_csv, sep=",", header=0)
    ids = list(df[column])
    if ids in ids:
        condition = df[column] == ids
        ind = str(df.index[condition])
        ind = "".join(re.findall(r"\[(\d+)\]", ind))
        return int(ind)
    else:
        return None


def set_loc_position(start, end, feat):
    """
    Set a new location for the feature in a genbank file.
    :param start: new star location of the gene.
    :param end: new ende loaction of the gene.
    :param feat: featrure to change the location.
    :return: feature with the new location.
    """
    start_pos = SeqFeature.ExactPosition(start)
    end_pos = SeqFeature.ExactPosition(end)
    my_location = SeqFeature.FeatureLocation(start_pos, end_pos, feat.location.strand)
    feat.location = my_location
    return feat.location


def set_compound_position(positions, feature):
    """
    Set a new loaction for a gen in a genbank file when the gene as a composed location, described in genbank files as join(....).
    :param positions: list of tuples where each tuple contains the new start ande the new end location in int [(new_start1, new_end1), (new_start2, new_end2)].
    :param feature: feature to change the location.
    :return: feature with the changed location.
    """
    new_positions = []
    for pos in positions:
        start, end = pos
        new_positions.append(SeqFeature.FeatureLocation(start, end, feature.location.strand))
    join_position = SeqFeature.CompoundLocation(new_positions)
    print(join_position)
    feature.location =  join_position
    return join_position

import numpy as np

def set_positions_2_by_2(positions):
    """
    positions: list with the positions to be grouped 2 by 2
    """
    print("\n positions" , positions)
    positions_2_by_2 = np.split(positions, np.arange(2,len(positions),2))
    print("\n grouped", positions_2_by_2)
    for lis in positions_2_by_2:
        print(lis)
    positions_2_by_2 = [(int(lis[0]), int(lis[1])) for lis in positions_2_by_2] 

    return positions_2_by_2


def get_locus_tag(feat, type_f):
    """
    Allows to retrieve the protein locus tag from the CDS of a genbank file.
    :param feat:
    :param type_f:
    :return:
    """
    if type_f == "CDS" or type_f == "gene":
        for k, v in feat.qualifiers.items():
            if k == "locus_tag":
                return v[0]
    else:
        return None


def get_protein_id(feat, type_f):
    """
    Allows to retrieve the protein id from the CDS of a genbank file.
    :param feat:
    :param type_f:
    :return:
    """
    if type_f == "CDS":
        for k, v in feat.qualifiers.items():
            if k == "protein_id":
                return v[0]
    else:
        return None


def get_protein_seq(feat, type_f):
    """
    Allows to retrieve the protein sequence from the CDS of a genbank file.
    :param feat: feature of the genbank file.
    :param type_f: type of the given feature
    :return: teh sequence of the protein from that feature
    """
    if type_f == "CDS":
        for k, v in feat.qualifiers.items():
            if k == "translation":
                return v[0]
    else:
        return None

def get_gene_abr(feat, type_f):
    """
    Allows to retive the gene abrevious from a CDS in genbank file.
    :param feat: feature from the genbank file
    :param type_f: feature type of the given feature
    :return: gene abreviation of the given feature
    """
    if type_f == "CDS" or type_f == "gene":
        for k, v in feat.qualifiers.items():
            if k == "gene":
                return v[0]
    else:
        return None


def write_id_file(file_name, info, path=None):
    """writes a txt file , from the given infromation.
    This function also allows to chose the path to save the files
     
    Args:
        file_name (str): string containing the name for the result file.
        info (list): list containing all the information to be placed in the file.
        path (str, optional): path to the folder to save the output files. If None the files will be saved in the current directory of the script.
    """
    if path == None:
        with open(file_name + ".txt", "w") as file:
            for ids in info:
                file.write(ids + "\n")
        file.close()
    else:
        with open(os.path.join(path, file_name + ".txt"), "w") as file:
            for ids in info:
                file.write(ids + "\n")
        file.close()


def write_file_txt_or_fasta(file_name, info, file_extension, path=None):
    """ Writes a txt file or a fasta file, dependeing on the chosen extension, from the given infromation.
    This function also allows to chose the path to save the files
     
    Args:
        file_name (str): string containing the name for the result file.
        info (list): list containing all the information to be placed in the file.
        file_extension (str): ".txt" or ".fasta" to choose which type of file will be generated.
        path (str, optional): path to the folder to save the output files. If None the files will be saved in the current directory of the script.
    """

    if path == None:
        with open(file_name + file_extension, "w") as file:
            for ids in info:
                file.write(ids + "\n")
        file.close()
    else:
        with open(os.path.join(path, file_name + file_extension), "w") as file:
            for ids in info:
                file.write(ids + "\n")
        file.close()    


def write_genes_fasta(header, seq, file_name):
    """ Writes a fasta file given an header, the sequence and the file name.

    Args:
        header (str): sting containing the wanted header to identify the sequence in the fasta
        seq (str): sequence to be placed in the fasta file
        file_name (str): name to the generated fasta file.
    """
    extension = ".fasta"
    file = file_name + extension
    with open(file, "a") as f:
        line = [header + "\n" + seq + "\n*2"]
        f.writelines(line)
    f.close()

