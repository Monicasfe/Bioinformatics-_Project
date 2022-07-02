from Bio import Entrez
from Bio.Seq import Seq
from Defs_Auxiliares import *
from All_csv_func import *
from Bio import SeqIO
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'

def get_proteins_ids(txt_files, searching_protein, add_to_file_name=None):
    genomes_dict = {}
    for i in range(len(txt_files)):
        with open(txt_files[i], "r") as f:
            txt_items = []
            for line in f:
                txt_items.append(line.split("\n")[0])
            genomes_dict[str(txt_files[i]).split(".")[0]] = txt_items

    Entrez.email = "example@mail.com"
    for k, v in genomes_dict.items():
        print(k)
        proteins_ids = []
        i = 0
        for v_in in v:
            genome_id = str(v_in).split(",")[-1].replace(" ", "")
            print(genome_id)
            genome = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb")
            seqs = open("Genomes.gb", "w")
            seqs.write(genome.read())
            phage = SeqIO.read("Genomes.gb", "gb")
            seqs.close()
            for feature in phage.features:
                if feature.type == "CDS":
                    for k1, v1 in feature.qualifiers.items():
                        if k1 == "product":
                            if v1[0] in searching_protein:
                                id_protein = get_protein_id(feature, feature.type)
                                if id_protein != None:
                                    proteins_ids.append(f"{v_in}, {str(id_protein)}")

            print("One is already done")

        print(k)

        if add_to_file_name != None:
            write_id_file(f"{add_to_file_name}{k}", proteins_ids)
        else: 
            write_id_file(f"{k}", proteins_ids)



# txt_file_names = ["NO_tail_protein.txt", "NO_baseplate_protein.txt", "NO_tail_fiber_protein.txt", "NO_tail_sheath_protein.txt"]

# get_proteins_ids(txt_files=txt_file_names, add_to_file_name="hypothetical_ports_ids_", searching_protein=["hypotetical protein"])

txt_files = ["NEW_tail_protein.txt", "NEW_baseplate_protein.txt", "NEW_tail_fiber_protein.txt", "NEW_tail_sheath_protein.txt"]

prots_names = [["tail", "tail protein"], ["baseplate", "baseplate protein"], ["tail fiber", "tail fiber protein"], ["tail sheath", "tail sheath protein"]]



def get_proteins_ids__222(txt_files, searching_protein, seach_from_list=False, add_to_file_name=None):
    """ txt_files list and seaching protein list must have have the order of the proteins
        
        EXEMPLE: txt_files = ["protein1.txt", "protein2.txt", "protein3.txt"]
                searching_protein = ["protein1", "protein2", "protein3"]

    Args:
        txt_files (_type_): _description_
        searching_protein (_type_): _description_
        seach_from_list (bool, optional): _description_. Defaults to False.
        add_to_file_name (_type_, optional): _description_. Defaults to None.
    """
    genomes_dict = {}
    for i in range(len(txt_files)):
        with open(txt_files[i], "r") as f:
            txt_items = []
            for line in f:
                txt_items.append(line.split("\n")[0])
            genomes_dict[str(txt_files[i]).split(".")[0].split("_")[-2]] = txt_items
    print(genomes_dict.keys())

    Entrez.email = "example@mail.com"

    i = 0

    for k, v in genomes_dict.items():
        print(k)
        proteins_ids = []
        
        for v_in in v:
            genome_id = str(v_in).split(",")[-1].replace(" ", "")
            genome_name = str(v_in).split(",")[0]
            print(genome_id)
            genome = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb")
            seqs = open("Genomes.gb", "w")
            seqs.write(genome.read())
            phage = SeqIO.read("Genomes.gb", "gb")
            seqs.close()
            duplicated_prots = {}
            for feature in phage.features:
                print(i)
                if feature.type == "CDS":
                    for k1, v1 in feature.qualifiers.items():
                        if k1 == "product":
                            if seach_from_list == True:

                                if any(k in pro for pro in searching_protein[i]):
                                    if v1[0] in searching_protein[i]:
                                        print("------------------------------------", v1[0])
                                        seq_protein = get_protein_seq(feature, feature.type)
                                        id_protein = get_protein_id(feature, feature.type)
                                        if id_protein and seq_protein!= None:
                                            duplicated_prots[id_protein] = len(str(seq_protein))

                            else:
                                if v1[0] in searching_protein:
                                    seq_protein = get_protein_seq(feature, feature.type)
                                    id_protein = get_protein_id(feature, feature.type)
                                    if id_protein and seq_protein!= None:
                                        duplicated_prots[id_protein] = len(str(seq_protein))

        

            print(duplicated_prots)
            if len(duplicated_prots.keys()) > 1:
                new_id_protein = max(duplicated_prots, key=duplicated_prots.get)
                print(new_id_protein)
                proteins_ids.append(f"{genome_id}, {genome_name}, {str(new_id_protein)}")
            else:
                proteins_ids.append(f"{genome_id}, {genome_name}, {str(list(duplicated_prots.keys())[0])}")

            print("One is already done")
        
        print(k)

        if add_to_file_name != None:
            write_id_file(f"{add_to_file_name}{k}",proteins_ids)
        else: 
            write_id_file(f"{k}",proteins_ids)

        i += 1

get_proteins_ids__222(txt_files=txt_files, searching_protein=prots_names, seach_from_list=True, add_to_file_name="NEW_NEW_prots_ids_")


def get_proteins_sequences(txt_files, add_to_file_name=None):
    genomes_dict = {}
    for i in range(len(txt_files)):
        with open(txt_files[i], "r") as f:
            txt_items = []
            for line in f:
                txt_items.append(line.split("\n")[0])
            genomes_dict[str(txt_files[i]).split(".")[0]] = txt_items

    Entrez.email = "example@mail.com"
    for k, v in genomes_dict.items():
        print(k)
        proteins_sequences = []
        i = 0
        for v_in in v:
            genome_id = str(v_in).split(",")[-1].replace(" ", "")
            print(genome_id)
            genome = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb")
            seqs = open("Genomes.gb", "w")
            seqs.write(genome.read())
            phage = SeqIO.read("Genomes.gb", "gb")
            seqs.close()
            for feature in phage.features:
                if feature.type == "CDS":
                    for k1, v1 in feature.qualifiers.items():
                        if k1 == "product":
                            if v1[0] == "hypothetical protein":
                                protein_id = get_protein_id(feature, feature.type)
                                sequence_protein = get_protein_seq(feature, feature.type)
                                print(protein_id)
                                if protein_id and sequence_protein != None:
                                    proteins_sequences.append(f">{protein_id}, {v_in}\n{sequence_protein}")
                            else:
                                print("Not an hypothetical protein")
        if add_to_file_name != None:
            write_file_txt_or_fasta(f"{add_to_file_name}{k}", proteins_sequences, file_extension=".fasta")
        else: 
            write_file_txt_or_fasta(f"{k}", proteins_sequences, file_extension=".fasta")          

# get_proteins_sequences(txt_files=txt_file_names, add_to_file_name="hypothetical_ports_sequences_")


def one_by_one_get_proteins_sequences(txt_files, add_to_file_name=None):
    genomes_dict = {}
    for i in range(len(txt_files)):
        with open(txt_files[i], "r") as f:
            txt_items = []
            for line in f:
                txt_items.append(line.split("\n")[0])
            genomes_dict[str(txt_files[i]).split(".")[0]] = txt_items

    Entrez.email = "example@mail.com"
    for k, v in genomes_dict.items():
        print(k)
        i = 0
        for v_in in v:
            genome_id = str(v_in).split(",")[-1].replace(" ", "")
            print(genome_id)
            proteins_sequences = []
            genome = Entrez.efetch(db="nucleotide", id=genome_id, rettype="gb")
            seqs = open("Genomes.gb", "w")
            seqs.write(genome.read())
            phage = SeqIO.read("Genomes.gb", "gb")
            seqs.close()
            for feature in phage.features:
                if feature.type == "CDS":
                    for k1, v1 in feature.qualifiers.items():
                        if k1 == "product":
                            if v1[0] == "hypothetical protein":
                                protein_id = get_protein_id(feature, feature.type)
                                sequence_protein = get_protein_seq(feature, feature.type)
                                print(protein_id)
                                if protein_id and sequence_protein != None:
                                    proteins_sequences.append(f">{protein_id}, {v_in}\n{sequence_protein}")
                                else:
                                    print("Not an hypothetical protein")
        
            if add_to_file_name != None:
                write_file_txt_or_fasta(f"{add_to_file_name}{k}_{v_in}", proteins_sequences, file_extension=".fasta", path="/home/acinom/Desktop/Mestrado/Projeto_Bioinf/Bioinf_project/genome_by_genome_sequences")
            else: 
                write_file_txt_or_fasta(f"{k}_{v_in}", proteins_sequences, file_extension=".fasta", path="/genome_by_genome_sequences")  
                            

# one_by_one_get_proteins_sequences(txt_files=txt_file_names, add_to_file_name="one_hypothetical_ports_sequences_")

