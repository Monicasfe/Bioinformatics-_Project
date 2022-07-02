import csv
from csv import writer
import pandas as pd
import os
import time


def creat_dataset(file_name, headers=True):
    """
    Creates a .csv file with column names
    :param file_name: String with the name that the file takes
    :param headers: List with column names
    :return: .cvs file empty
    """
    for i in file_name:
        if headers == False:
            with open(i + ".csv", "w", encoding="UTF8", newline="") as file:
                writer = csv.writer(file, delimiter=",")
        else:
            with open(i + ".csv", "w", encoding= "UTF8", newline="") as file:
                writer = csv.writer(file, delimiter = ",")
                writer.writerow(headers)

def append_row(file, list):
    """
    Appends a row to csv file.
    :param file: file to add the row
    :param list: list with the information
    :return:
    """
    with open(file, 'a+', encoding= "UTF8", newline='') as dataset:
        csv_writer = writer(dataset)
        csv_writer.writerow(list)


def join_csv(file, new_file, extension, remove=False, indexes=None):
    """
    Joins to csv files.
    :param file: list of files to concatenate.
    :param new_file: str with the name of the output file
    :param extension: ".csv"
    :param remove: bolean, if True remove the indexes in the indexes list. By default is False.
    :param indexes: list containing the indexes to remove if "remove=True". By default is None
    :return:
    """
    all_in_one = []
    nfile_name = new_file + extension
    for f in file:
        print(f)
        fi = pd.read_csv(str(f))
        print(fi)
        all_in_one.append(fi)
    concat = pd.concat(all_in_one, ignore_index=True)
    if remove == True:
        for i in indexes:
            del concat[i]
        final_concat = concat.to_csv(nfile_name, index=False, encoding = "utf-8-sig")
        df = pd.read_csv(nfile_name, header=0)
        print(df)
        return final_concat
    else:
        final_concat_2 = concat.to_csv(nfile_name, index=False, encoding="utf-8-sig")
        df = pd.read_csv(nfile_name, header=0)
        print(df)
        return final_concat_2

# data = ["TerL_Genomes", "No_TerL_Genomes", "Caquinha", "Caquinha_2"]
#
# headers = ["Accession_Number", "Phage_name", "Host", "Locus", "Start_protein", "Size", "DNA"]
#
# start_time = time.time()

headers_2 = ["Accession_Number", "Phage_name", "Host", "Size", "DNA"]
#
# creat_dataset(["Baseplate_All_phages"], headers_2)
#
# if __name__ == "__main__":
#     print("--- %s seconds ---" % (time.time() - start_time))