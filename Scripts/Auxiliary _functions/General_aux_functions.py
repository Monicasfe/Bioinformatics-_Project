from hashlib import new
import re
from nbformat import write
import pandas as pd
from Bio import SeqFeature
import os


def remove_id(list, ids):
    remove_id = []
    for i in list:
        if ids in i:
            remove_id.append(i)
    for j in remove_id:
        if j in list:
            list.remove(j)
    return list



def sort_by_id(list, ids):
    new = []
    for i in list:
        if ids in i:
            new.insert(0, i)
        if ids not in i:
            new.append(i)
    return new


def get_feature(seqio_obje, feature):
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
    start_pos = SeqFeature.ExactPosition(start)
    end_pos = SeqFeature.ExactPosition(end)
    my_location = SeqFeature.FeatureLocation(start_pos, end_pos, feat.location.strand)
    feat.location = my_location
    return feat.location


def set_compound_position(positions, feature):
    new_positions = []
    for pos in positions:
        start, end =  pos
        new_positions.append(SeqFeature.FeatureLocation(start, end, feature.location.strand))
    join_position = SeqFeature.CompoundLocation(new_positions)
    print(join_position)
    feature.location =  join_position
    return join_position

import numpy as np

def set_positions_2_by_2(positions):
    """
    positin: list with the positions to be grouped 2 by 2
    """
    print("\n fuking positions" , positions)
    positions_2_by_2 = np.split(positions, np.arange(2,len(positions),2))
    print("\n merda agrupada", positions_2_by_2)
    for lis in positions_2_by_2:
        print(lis)
    positions_2_by_2 = [(int(lis[0]), int(lis[1])) for lis in positions_2_by_2] 

    return positions_2_by_2


def get_locus_tag(feat, type_f):
    if type_f == "CDS" or type_f == "gene":
        for k, v in feat.qualifiers.items():
            if k == "locus_tag":
                return v[0]
    else:
        return None


def get_protein_id(feat, type_f):
    if type_f == "CDS":
        for k, v in feat.qualifiers.items():
            if k == "protein_id":
                return v[0]
    else:
        return None


def get_protein_seq(feat, type_f):
    if type_f == "CDS":
        for k, v in feat.qualifiers.items():
            if k == "translation":
                return v[0]
    else:
        return None

def get_gene_abr(feat, type_f):
    if type_f == "CDS" or type_f == "gene":
        for k, v in feat.qualifiers.items():
            if k == "gene":
                return v[0]
    else:
        return None


def write_id_file(file_name, info, path=None):
    """_summary_

    Args:
        file_name (str): _description_
        info (list): _description_
        path (_type_, optional): _description_. Defaults to None.
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
    extension = ".fasta"
    file = file_name + extension
    with open(file, "a") as f:
        line = [header + "\n" + seq + "\n*2"]
        f.writelines(line)
    f.close()

