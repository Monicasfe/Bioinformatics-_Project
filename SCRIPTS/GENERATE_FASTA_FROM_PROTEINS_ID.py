from importlib.resources import path
from Bio import Entrez
from Bio.Seq import Seq
from more_itertools import only
from SCRIPTS.AUXLILIARY_FUNCTIONS.GENERAL_AUX_FUNCTIONS import *
from SCRIPTS.AUXLILIARY_FUNCTIONS.AUX_CSV_FUNCTIONS import *
from Bio import SeqIO
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'


def get_fasta_from_prot_id(txt_files, add_to_file_name=None, all_in_one_file=True, path_files=None):
    """
    This function will generate a .fasta file containing all the genes and their position extracted form a genbank (gb) file, given an Accesion Number and a protein id.    
    Each given gb file will give origin to .cvs file.
    This function accepts only one file or a directory containing more than one gb file
    :param dir: string of the directory to the gb file/s
    :return: one csv file per gb file with the name of the organism from the gb file
    """

    for j in range(len(txt_files)):
        txt_dict = {}

        file_name = str(txt_files[j]).split(".")[0]

        with open(txt_files[j], "r") as f:
            for line1 in f:
                txt_item = " ".join(line1.split(", ")[1:]).replace("\n", "").replace(" ",", ")
                txt_dict[str(line1.split(",")[0]).split(".")[0]] = txt_item
            
            Entrez.email = "example@mail.com"
            prots_sequences = []
            for k, v in txt_dict.items():

                only_one_seq = []

                print(k)

                genome = Entrez.efetch(db="nucleotide", id=k, rettype="gb")
                seqs = open("Genomes.gb", "w")
                seqs.write(genome.read())
                seqs.close()
                phage = SeqIO.read("Genomes.gb", "gb") 

                sequence_genome = Seq(phage.seq).lower()
                
                print(len(sequence_genome))

                protein_id = str(v).split(", ")[-1]
                print(protein_id)
                for feature in phage.features:
                    if feature.type == "CDS":
                        for k1, v1 in feature.qualifiers.items():
                            if k1 == "product":
                                prot_id = get_protein_id(feature, feature.type)
                                if protein_id == prot_id:
                                    start_loc = int(feature.location.start) #-1 por causa do zero
                                    print("start_loc", start_loc)
                                    end_loc = int(feature.location.end)
                                    print("end_loc", end_loc)
                                    gene_loc = (start_loc, end_loc)
                                    strand = int(str(feature.location.strand))
                                    print(strand)
                                    if strand == 1:
                                        seq = sequence_genome[start_loc:end_loc]
                                    else:
                                        seq = sequence_genome[start_loc:end_loc].complement() 


                                    if all_in_one_file==True:
                                        prots_sequences.append(f">{k}, {v}, {v1[0]}, {gene_loc} \n{seq} \n")
    
                                        print(prots_sequences)

                                        print("\n")

                                    else:
                                        only_one_seq.append(f">{k}, {v}, {v1[0]}, {gene_loc} \n{seq} \n")

                                        if add_to_file_name != None:
                                            if path_files!=None:
                                                write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}_{k}", info=only_one_seq, file_extension=".fasta", path=path_files[j])
                                            else:
                                                write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}_{k}", info=only_one_seq, file_extension=".fasta")
                                        else: 
                                            if path_files!=None:
                                                write_file_txt_or_fasta(file_name=f"{file_name}_{k}", info=only_one_seq, file_extension=".fasta", path=path_files[j])
                                            else:
                                                write_file_txt_or_fasta(file_name=f"{file_name}_{k}", info=only_one_seq, file_extension=".fasta")

        
        

                print("\n")

        if all_in_one_file==True:
            if add_to_file_name != None:
                if path_files!=None:
                    write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=prots_sequences, file_extension=".fasta", path=path_files[j])
                else:
                    write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=prots_sequences, file_extension=".fasta")
            else: 
                if path_files!=None:
                    write_file_txt_or_fasta(file_name=f"{file_name}", info=prots_sequences, file_extension=".fasta", path=path_files[j])
                else:
                    write_file_txt_or_fasta(file_name=f"{file_name}", info=prots_sequences, file_extension=".fasta", path=path_files[j])
        else:
            print("One is Done")
    


# prots_files = ["NEW_NEW_prots_ids_baseplate.txt", "NEW_NEW_prots_ids_tail.txt", "NEW_NEW_prots_ids_fiber.txt", "NEW_NEW_prots_ids_sheath.txt"]
#
# paths = ["nt_fasta_sequences/Baseplate", "nt_fasta_sequences/Tail", "nt_fasta_sequences/Tail_fiber", "nt_fasta_sequences/Tail_sheath"]
#
# get_fasta_from_prot_id(txt_files=prots_files, add_to_file_name="FASTA_", all_in_one_file=False, path_files=paths)


def get_fasta_from_prot_location(csv_files, add_to_file_name=None, all_in_one_file=True, path_files=None):
    """
    This function will generate a .fasta file containing all the genes and their position extracted form a genbank (gb) file, given an Accesion Number and a protein id.    
    Each given gb file will give origin to .cvs file.
    This function accepts only one file or a directory containing more than one gb file
    :param dir: string of the directory to the gb file/s
    :return: one csv file per gb file with the name of the organism from the gb file
    """

    for i in range(len(csv_files)):

        file_name = "_".join(str(csv_files[i]).split(".")[0].split("_")[2:])

        print(file_name)

        prots_sequences = []
        
        df = pd.read_csv(csv_files[i], sep=",", header=0)
        print(df.shape[0])
        for j in range(df.shape[0]):
            only_one_seq = []
            ids = df["Accession"][j]
            start_loc = int(str(df["location"][j]).split(",")[0])
            end_loc = int(str(df["location"][j]).split(",")[1])
            name = df["Common_Name"][j]
            gene_loc = (start_loc, end_loc)

            Entrez.email = "example@mail.com"
           

            print(ids)

            genome = Entrez.efetch(db="nucleotide", id=ids, rettype="gb")
            seqs = open("Genomes.gb", "w")
            seqs.write(genome.read())
            seqs.close()
            phage = SeqIO.read("Genomes.gb", "gb") 

            sequence_genome = Seq(phage.seq).lower()
            
            print(len(sequence_genome))

            seq = sequence_genome[start_loc:end_loc]
            
            if all_in_one_file==True:
                prots_sequences.append(f">{ids}, {name}, {file_name}, {gene_loc} \n{seq} \n")
    
                print(prots_sequences)

                print("\n")

            else:
                only_one_seq.append(f">{ids}, {name}, {file_name}, {gene_loc} \n{seq} \n")

                if add_to_file_name != None:
                    if path_files!=None:
                        write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}_{ids}", info=only_one_seq, file_extension=".fasta", path=path_files[i])
                    else:
                        write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}_{ids}", info=only_one_seq, file_extension=".fasta")
                        
                else: 
                    if path_files!=None:
                        write_file_txt_or_fasta(file_name=f"{file_name}_{ids}", info=only_one_seq, file_extension=".fasta", path=path_files[i])
                    else:
                        write_file_txt_or_fasta(file_name=f"{file_name}_{ids}", info=only_one_seq, file_extension=".fasta")
                        
        
        if all_in_one_file==True:
            if add_to_file_name != None:
                if path_files!=None:
                    write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=prots_sequences, file_extension=".fasta", path=path_files[i])
                else:
                    write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=prots_sequences, file_extension=".fasta")
            else: 
                if path_files!=None:
                    write_file_txt_or_fasta(file_name=f"{file_name}", info=prots_sequences, file_extension=".fasta", path=path_files[i])

                else:
                    write_file_txt_or_fasta(file_name=f"{file_name}", info=prots_sequences, file_extension=".fasta")
        else:
            print("One is done")

# csv_files = ["selected_homology_baseplate_protein.csv", "selected_homology_tail_protein.csv",  "selected_homology_tail_fiber_protein.csv", "selected_homology_tail_sheath_protein.csv"]
#
# get_fasta_from_prot_location(csv_files, all_in_one_file=False, path_files=paths)
