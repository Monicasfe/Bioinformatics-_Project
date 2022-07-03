from Bio import Entrez
from Bio.Seq import Seq
from SCRIPTS.AUXLILIARY_FUNCTIONS.GENERAL_AUX_FUNCTIONS import *
from SCRIPTS.AUXLILIARY_FUNCTIONS.AUX_CSV_FUNCTIONS import *
import time
import os
from Bio import SeqIO
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

def check_proteins_in_dataset(dataset, proteins, phages_file, prots_counts, add_to_file_name=None):
    """ This funtions' output is a file with all the genomes that contain the wanted protein notated in their genbank file. 

    Args:
        dataset (csv): A cvs file containing the accession number and the genomes' name. 
        proteins (list): Alist containing all the version of the names for the wanted protein to search.
        phages_file (str): Str by which the output file will be named of.
        prots_counts (int): Int that works as a counter to check if the genome has all the wanted proteins. 
    """
    Entrez.email = "example@mail.com"
    df = pd.read_csv(dataset, sep=",", header=0)
    ids = list(df["Accession_Number"])
    names = list(df["Phage_name"])
    have_all_prots = []
    i=0
    print(len(ids))
    while i < len(ids):
        prots_names = []
        n_prots = 0
        print(n_prots)
        j = ids[i]
        print(j, i)
        genome = Entrez.efetch(db="nucleotide", id=j, rettype="gb")
        seqs = open("Genomes.gb", "w")
        seqs.write(genome.read())
        phage = SeqIO.read("Genomes.gb", "gb")
        seqs.close()
        for feature in phage.features:
            j = 0
            if feature.type == "CDS":
                for k, v in feature.qualifiers.items():
                    if k == "product":
                        if v[0] in proteins:
                            if prots_names == []:
                                    print(v[0])
                                    n_prots += 1
                                    print(n_prots)
                                    prots_names.append(str(v[0]))
                            else:
                                # if "putative" in str(v[0]).split():
                                print("------------" + str(v[0]))
                                for name in prots_names:
                                    if " ".join(str(v[0]).split()[-2:]) in name:
                                        print(" ".join(str(v[0]).split()[-2:]))
                                        j+=1
                                if j == 1:
                                    print(f"{v[0]} já existe")

                                else:
                                    print(v[0])
                                    n_prots += 1
                                    print(n_prots)
                                    prots_names.append(str(v[0]))

        if n_prots == prots_counts:
            have_all_prots.append(f"{names[i]}, {ids[i]}")
            print("aqui")
        i += 1
        # if i == len(ids) +1:
        #     print(i)
        #     break
    if have_all_prots != []:
        if add_to_file_name != None:
            write_id_file(f"{add_to_file_name}{phages_file}", have_all_prots)
    else:
        print("list is empty")


# baseplate = ["baseplate protein", "baseplate"]

# tail_sheath = ["tail sheath protein", "tail sheath"]

# tail_fiber = ["tail fiber protein", "tail fiber"]

# tail = ["tail protein"]


# all_searches = [baseplate, tail, tail_fiber, tail_sheath]


# for i in all_searches:
    # if len(i) > 3:
        # check_proteins_in_dataset("NEW_All_phages.csv", i, "All_proteins", prots_counts=4)
    # else:
        # check_proteins_in_dataset("NEW_All_phages.csv", i, str(i[0]).replace(" ", "_"), prots_counts=1, add_to_file_name="NEW_NEW_")

            
