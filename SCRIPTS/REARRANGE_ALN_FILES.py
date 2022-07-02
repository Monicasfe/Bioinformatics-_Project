import re
import numpy as np
from Defs_Auxiliares import *

def arrange_aln_files_to_fasta_style(aln_files, add_to_file_name=None, path_files=None):
    
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
            # print(lines)

            # lines = lines[3:]
            # print(len(lines))
            # print(lines[-2])
            # for j in range(len(lines)):
            #     if len(lines) == 0:
            #         break
            #     else:
            #         print(len(lines))
            #         print(lines)
            #         j = j-j
            #         seq = []
            #         # print(lines[j])
            #         splited = lines[0].split(" ")
            #         id_number = splited[0]
            #         seq.append(splited[-1])
            #         lines.remove(lines[j])
            #         print(id_number)
                    
            #         # for line_3 in lines:
            #         #     if id_number in line_3.split(" ")[0]:
                            
            #         # if any(id_number in x for x in lines):
            #         idx = [ids for ids, l in enumerate(lines) if id_number in l.split(" ")[0]]
            #         print(idx)
            #         for h in idx:
            #             print(lines[h].split(" ")[0])
            #             seq.append(lines[h].split(" ")[-1])
                    
            #         for g in idx:
            #             # print(len(lines))
            #             lines.remove(lines[int(g)])
                    
            #         # print(seq)
            #         aln_dict[id_number] = seq

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

        sequences_aln = []

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