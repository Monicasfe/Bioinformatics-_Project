from Bio import Entrez
from Bio import SeqIO
from AUXLILIARY_FUNCTIONS.GENERAL_AUX_FUNCTIONS import *
import http.client
http.client.HTTPConnection._http_vsn = 10
http.client.HTTPConnection._http_vsn_str = 'HTTP/1.0'
import urllib.request
from urllib.error import HTTPError

def generate_dataset(search_query, db_s, file_1):
    """
    Generates a datasets containing all the organisms found the NCBI database from a given searching query.

    This fucntion has the requierement of a specific header for the dataset whihc is:

    "Accession_number", "Organism name", "Host", "Sequence_length"

    Inside this function there are filter terms that can be changed according to users necessity.
    :param search_query: String containing the search query to provide to the database.
    :param db_s: NCBI database, nucleotide, protein and so on.
    :param file_1: CSV file where all the good genomes and their information will be written.
    :return: To csv files containing genomes information
    """
    Entrez.email = "example@mail.com"
    handle = Entrez.esearch(db=db_s, retmax=800, term=search_query, idtype="acc")
    record = Entrez.read(handle)
    ids = []
    dataset = []
    for id in record["IdList"]:
        ids.append(id)
    print(ids)
    ids = remove_id(ids, "NZ_CP")
    print(ids)
    print(len(ids))
    i = 0
    phages_names = []
    while i < len(ids):
        print("i initial", i)
        j = ids[i]
        genome = Entrez.efetch(db="nucleotide", id=j, rettype="gb")
        seqs = open("Genomes.gb", "w")
        seqs.write(genome.read())
        phage = SeqIO.read("Genomes.gb", "gb")
        if get_feature(phage, "CDS") is None:
            i += 1
            continue
        description = phage.description
        if "UNVERIFIED" in description:
            i += 1
            continue
        elif ("phage" or "prophage") not in description:
            i += 1
            continue
        elif "Staphylococcus" not in description:
            i += 1
            continue
        elif "partial genome" in description:
            i += 1
            continue
        elif "host-range" in description:
            i += 1
            continue
        else:
            print(description)
        source = get_feature(phage, "source")
        CDS = get_feature(phage, "CDS")
        if CDS is None:
            i += 1
            continue
        for k, v in source.qualifiers.items():
            if k == "organism":
                p_name = v[0]
                if ("phage" or "virus") not in p_name:
                    i += 1
                    continue
                else:
                    print(p_name)
                host = p_name.split()[0]
            if k == "host":
                host_1 = v[0]
                host = host_1
                print(host)
        if "Staphylococcus" in p_name:
            p_name = " ".join(p_name.split()[2:])#[0] tirei porque quero o que há depois do 1
            print(p_name)
        else:
            print(p_name)
        if p_name not in phages_names:
            phages_names.append(p_name)
        else:
            print("!!!!!!!Duplicated removed!!!!!!")
            i += 1
            continue
        print(phages_names)
        id = phage.name
        length = len(phage.seq)
        # dataset = [id, p_name, host, length, seq]
        # append_row(file_1, dataset)
        dataset.append(f"{id},{p_name},{host}, {int(length)}")
        print(dataset)
        i += 1
        print("number of phages", len(phages_names))
        print(i)

    write_id_file(file_name=file_1, info=dataset)

# searching_query = "staphylococcus phage[All Fields] AND complete genome[All Fields]"

# file_1 = "NEW_All_phages.csv"
# file_1 = "NEW_All_phages"

# start_time = time.time()

# generate_dataset(searching_queryquery, "nucleotide", file_1)

# if __name__ == "__main__":
    # print("--- %s seconds ---" % (time.time() - start_time))
