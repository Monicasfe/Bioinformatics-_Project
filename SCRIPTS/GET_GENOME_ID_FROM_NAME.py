from Bio import Entrez
from Bio.Seq import Seq
from numpy import dtype
from SCRIPTS.AUXLILIARY_FUNCTIONS.GENERAL_AUX_FUNCTIONS import *
from SCRIPTS.AUXLILIARY_FUNCTIONS.AUX_CSV_FUNCTIONS import *
from Bio import SeqIO
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

def get_genome_id_from_name(dataset, csv_name_files, add_to_file_name):
    """
    From a given dataset having the accession number and name of all the genomes we are working with it gets teh accession number of the genome based on their name.
    :param dataset: initial dataset in a csv format, having among other inormation teh accession number and genome name.
    :param csv_name_files: csv files having the name of the phages from which we want their accession number.
    :param add_to_file_name: str having a prefix to add to the name of teh output file.
    :return:
    """
    df = pd.read_csv(dataset, sep=",", header=0)

    print(df.shape)

    ids_names_dict ={}

    for i in range(df.shape[0]):
        ids_names_dict[df["Accession_Number"][i]] = df["Phage_name"][i]

    print(ids_names_dict)

    for j in range(len(csv_name_files)):
        file_name = " ".join(str(csv_name_files[j]).split(".")[0].split("_")[-2:]).replace(" ", "_")
        lines_info = []
        data_df = pd.read_csv(csv_name_files[j], sep=",", header=0)
        print(data_df)
        for l in range(len(data_df)):
            name = data_df["Common_Name"][l]
            protein = data_df["Accession"][l]
            if name in ids_names_dict.values():
                # print(data_df["Common_Name"][l])     
                key = [k for k, v in ids_names_dict.items() if v == data_df["Common_Name"][l]]
                lines_info.append(f"{key[0]}, {name}, {protein}")
            else:
                print(f"Does not exist {name}")

        if add_to_file_name != None:
            write_id_file(f"{add_to_file_name}{file_name}", lines_info)
        else: 
            write_id_file(f"{file_name}", lines_info)




name_files = ["selected_homology_baseplate_proteins.csv", "selected_homology_tail_fiber_proteins.csv", "selected_homology_tail_proteins.csv", "selected_homology_tail_sheath_proteins.csv"]


get_genome_id_from_name("NEW_All_phages.csv", csv_name_files=name_files, add_to_file_name="id_name_homology_")