from audioop import add
from Defs_Auxiliares import *
from All_csv_func import *
import os
import pandas as pd


def select_genomes_no_prots(dataset, txt_files, add_to_file_name=None):
    df = pd.read_csv(dataset, sep=",", header=0)
    ids = list(df["Accession_Number"])
    proteins_dict = {}
    for i in range(len(txt_files)):
        with open(txt_files[i], "r") as f:
            txt_items = []
            for line in f:
                txt_items.append(line.split(", ")[1].replace("\n", ""))
            proteins_dict[str(txt_files[i]).split(".")[0]] = txt_items
    print(proteins_dict)

    no_notation_genomes_dict = {}

    for k, v in proteins_dict.items():
        no_notation_genomes = []
        print(v)
        for id_genome in ids:
            print(id_genome)
            if id_genome in v:
                print("This one exists")
                
            else: 
                no_notation_genomes.append(id_genome)
                print("A new one")
        no_notation_genomes_dict[k] = no_notation_genomes

    for k1, v1 in no_notation_genomes_dict.items():
        if add_to_file_name is None:
            write_id_file(k1, v1) 
        else: 
            write_id_file(file_name=f"{add_to_file_name}{k1}", info=v1)   


txt_list = ["NEW_tail.txt", "NEW_baseplate_protein.txt", "NEW_tail_fiber_protein.txt", "NEW_tail_sheath_protein.txt"]
# out_file_names = ["NO_tail_protein.txt", "NO_baseplate_protein.txt", "NO_tail_fiber_protein.txt", "NO_tail_sheath_protein.txt"]


select_genomes_no_prots("NEW_All_phages.csv", txt_list, add_to_file_name="NO_")

            