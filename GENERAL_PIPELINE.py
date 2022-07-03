from SCRIPTS.AUXLILIARY_FUNCTIONS.AUX_CSV_FUNCTIONS import creat_dataset
from SCRIPTS.GENERATE_PHAGE_DATASET import *
from SCRIPTS.CHECK_PROTEINS import *
from SCRIPTS.GET_GENOMES_WITHOUT_THE_PROTEINS import *
from SCRIPTS.GET_PROTEINS_ID import *
from SCRIPTS.GENERATE_FASTA_FROM_PROTEINS_ID import *
from SCRIPTS.ALIGNMENTS import alignment
from SCRIPTS.REARRANGE_ALN_FILES import *

start_time = time.time()

###### CREATE a csv file to allocate the dataset information ###

headers = ["Accession_Number", "Phage_name", "Host", "Size"]

creat_dataset(["NEW_All_phages.csv"], headers)

################################### PART 1 #################################
    ## Generation of the dataset by giving a searchin query, a csv files and the NCBI database type ##

query = "staphylococcus phage[All Fields] AND complete genome[All Fields]"

file_1 = "NEW_All_phages.csv"

generate_dataset(query, "nucleotide", file_1)

################################### PART 2 #################################

 ##  Check which genomes have all the wanted proteins and which genomes have each protein. Generates a file per protein and one for all proteins combined ##

baseplate = ["baseplate protein", "baseplate"]

tail_sheath = ["tail sheath protein", "tail sheath"]

tail_fiber = ["tail fiber protein", "tail fiber"]

tail = ["tail protein"]


all_searches = [baseplate, tail, tail_fiber, tail_sheath]

for i in all_searches:
   if len(i) > 3:
       check_proteins_in_dataset("NEW_All_phages.csv", i, "All_proteins", prots_counts=4)
   else:
       check_proteins_in_dataset("NEW_All_phages.csv", i, str(i[0]).replace(" ", "_"), prots_counts=1, add_to_file_name="NEW")


################################### PART 3 #################################
    ## get from the original dataset wich genomes do not have teh wanted proteins. This is performed by each of teh proteins.

txt_list = ["NEW_tail_protein.txt", "NEW_baseplate_protein.txt", "NEW_tail_fiber_protein.txt", "NEW_tail_sheath_protein.txt"]

select_genomes_no_prots("NEW_All_phages.csv", txt_list, add_to_file_name="NO_")



################################### PART 4 #################################
    ## Get the ids of the wanted proteins inside the genome ##
        ## laong with this step you shoul use teh files with the genomes without the wanted proteins to run tblastn on NCBI to obatin homologous proteins ans extract the information of genome name, accession number and protein location inside the genome (start, end)##

txt_files = ["NEW_tail_protein.txt", "NEW_baseplate_protein.txt", "NEW_tail_fiber_protein.txt", "NEW_tail_sheath_protein.txt"]

prots_names = [["tail", "tail protein"], ["baseplate", "baseplate protein"], ["tail fiber", "tail fiber protein"], ["tail sheath", "tail sheath protein"]]

get_proteins_ids(txt_files=txt_files, searching_protein=prots_names, seach_from_list=True, add_to_file_name="NEW_prots_ids_")


################################### PART 5 #################################

    ## generate a fasta for each group of porteins ##
        ## using the proteins ids ##

prots_files = ["NEW_prots_ids_baseplate.txt", "NEW_prots_ids_tail.txt", "NEW_prots_ids_fiber.txt", "NEW_prots_ids_sheath.txt"]

paths = ["nt_fasta_sequences/Baseplate", "nt_fasta_sequences/Tail", "nt_fasta_sequences/Tail_fiber", "nt_fasta_sequences/Tail_sheath"]

get_fasta_from_prot_id(txt_files=prots_files, add_to_file_name="FASTA_", all_in_one_file=False, path_files=paths)

        ## using proteins from homology ##

csv_files = ["selected_homology_baseplate_protein.csv", "selected_homology_tail_protein.csv",  "selected_homology_tail_fiber_protein.csv", "selected_homology_tail_sheath_protein.csv"]

get_fasta_from_prot_location(csv_files, all_in_one_file=False, path_files=paths)

    ## After this step join both fasta files from each protein into one

################################### PART 6 #################################
    ## Perform the alignments using ClustalW ##


diretoria = r"/path/to/tool/clustalw2"

dir_file = [r"/path/to/file/Fasta_baseplate.fasta", r"/path/to/file/Fasta_tail.fasta", r"/path/to/file/Fasta_tail_fiber.fasta", r"/path/to/file/Fasta_tail_sheath.fasta"]
out_file = ["Baseplate_aln", "Tail_aln", "Tail_fiber_aln", "Tail_sheath_aln"]

for i in range(len(dir_file)):
    alignment(diretoria, dir_file[i], output_file=out_file[i])

################################### PART 7 #################################

    ## Conversion of the structure of the aln files into a fasta sctruture to allow the converion of teh aln file to a vcf file. Ti generate the vcf file it was used the SNP-sites tool

    aln_files = ["Fasta_baseplate.aln", "Fasta_tail.aln", "Fasta_tail_fiber.aln", "Fasta_tail_sheath.aln"]

    paths_files = ["SNP/baseplate", "SNP/tail", "SNP/tail_fiber", "SNP/tail_sheath"]

    arrange_aln_files_to_fasta_style(aln_files=aln_files, add_to_file_name="Like_fasta_", path_files=paths_files)

################################### PART 8 #################################
    ## To perform SNP caling simply run each script ##


################################### PART 9 #################################
    ## To draw the calsograms junt go over the jupyter notebook and run it. ##


if __name__ == "__main__":
    print("--- %s seconds ---" % (time.time() - start_time))