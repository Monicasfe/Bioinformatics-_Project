import re
import numpy as np
from SCRIPTS.AUXLILIARY_FUNCTIONS.GENERAL_AUX_FUNCTIONS import *

def arrange_aln_files_to_fasta_style(aln_files, add_to_file_name=None, path_files=None):
    """
    Converts an alignment file sln structure to on alignment file with the structure of a fasta file.
    :param aln_files: list with the aln files to convert.
    :param add_to_file_name: prefix to add to the output file name.
    :param path_files: path to where the output files should be saved. By default is None, meaning that the output file will be saved in te corrent folder
    :return:
    """
    
    for i in range(len(aln_files)):
        with open(aln_files[i], "r") as f:
            aln_dict = {}
            lines = []

            for line in f:
                if re.findall(r'([A-Z])', line)!=[]:
                    line_1 = [string for string in line.split(" ") if string != ""]
                    line_2 = " ".join(line_1)
                    lines.append(str(line_2).replace("\n", ("")))
                    # print(line_2)
                else:
                    print("esta merda existe", line)
            lines = lines[1:]

            ids = []
            for line in lines:
                    line = line.split(" ")  # dividir entre ID e SEQ
                    ids.append(line[0])
                    if line[0] not in aln_dict:
                        aln_dict[line[0]] = line[1]
                    else:
                        aln_dict[line[0]] = aln_dict[line[0]] + line[1]

        
        print(len(aln_dict), len(np.unique(ids)), len(aln_dict.values()))
        print(aln_dict.keys())

        # print(aln_dict)

        for k, v in aln_dict.items():
            k = k.replace("_", "")
            sequences_aln.append(f">{k} \n{v} \n")
        file_name = "_".join(aln_files[i].split(".")[:1])
        if add_to_file_name != None:
            if path_files != None:
                write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=sequences_aln, file_extension=".aln", path=path_files[i])
            else:
                write_file_txt_or_fasta(file_name=f"{add_to_file_name}{file_name}", info=sequences_aln, file_extension=".aln")
        else:
            if path_files != None:
                write_file_txt_or_fasta(file_name=f"{file_name}", info=sequences_aln, file_extension=".aln", path=path_files[i])
            else:
                write_file_txt_or_fasta(file_name=f"{file_name}", info=sequences_aln, file_extension=".aln")


aln_files = ["Fasta_baseplate.aln", "Fasta_tail.aln", "Fasta_tail_fiber.aln", "Fasta_tail_sheath.aln"]

paths_files = ["SNP/baseplate", "SNP/tail", "SNP/tail_fiber", "SNP/tail_sheath"]

arrange_aln_files_to_fasta_style(aln_files=aln_files, add_to_file_name="Like_fasta_", path_files=paths_files)