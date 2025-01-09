import pandas
import os
from pathlib import Path
import re

hla_directory = "~/genotyping_data/hisat_hla_results/"
kir_directory = "~/genotyping_data/unadjusted_kir_genotyping_results/"

#get a list of all the files in the directories
hla_file_list = []
kir_file_list = []

#add the files to a list
for file in os.scandir(hla_directory):
    hla_file_list.append(Path(file))
for file in os.scandir(kir_directory):
    kir_file_list.append(Path(file))

files = pandas.DataFrame(list(zip(hla_file_list, kir_file_list)), columns = ["HLA_files","KIR_files"])

#Check that each numbered file has it's correct pair

def fasta_files(file_name):
    sequences = []
    with open(file_name, "r") as file:
        current_sequence = ""
        for line in file:
            line = line.strip()  #remove newlines
            if line.startswith(">"):  #this is the sequence label
                if current_sequence:  #if there is a previous sequence, add it to the list
                    sequences.append(current_sequence)
                current_sequence = [line]  #start a new sequence
            else:
                current_sequence.append(line)  #add line to the current sequence
        sequences.append(current_sequence)  #add the last sequence after exiting the loop
    #combine sequences into single strings
    combined_sequences = []
    for sequence in sequences:
        combined_sequences.append([sequence[0], ''.join(sequence[1:])])

    return combined_sequences

#read in the files, put them into their lists, this is used later on
HLA_A_sequences = fasta_files("~/project2_KIR/A_prot.fasta")
HLA_B_sequences = fasta_files("~/project2_KIR/B_prot.fasta")
HLA_C_sequences = fasta_files("~s/project2_KIR/C_prot.fasta")

def extract_patient_alleles(hla_file, kir_file):
    #get HLA files into dataframe
    hla_filename = os.path.basename(str(hla_file))
    hla_sequencing_id = hla_filename.removeprefix("result_LAI550A").removesuffix("_hla.csv")
    hla_genotype_table = pandas.read_csv(hla_file)
    
    #get KIR files into dataframe
    kir_filename = os.path.basename(str(kir_file))
    kir_sequencing_id = kir_filename.removeprefix("kir_run_LAI550A").removesuffix("_genotype.tsv")
    kir_genotype_table = pandas.read_table(kir_file, index_col=0, names=["gene_name","num_alleles","allele_1","abundance_1","qscore_1","allele_2","abundance_2","qscore_2"])
    
    #this pairs the files
    if kir_sequencing_id == hla_sequencing_id:
        sequencing_id = kir_sequencing_id

        #get the HLA alleles
        pattern = r'\ \(abundance:\s\d+(\.\d+)?%\)'
        for idx, row in hla_genotype_table.iterrows():
        #A alleles
            HLA_A_1 = re.sub(pattern, "", hla_genotype_table.at[0,"EM: A"])
            if idx > 0: #make sure that index 1 isn't out of range
                if hla_genotype_table.at[1,"EM: A"] == "None":
                    HLA_A_2 = HLA_A_1
                else:
                    HLA_A_2 = re.sub(pattern, "", hla_genotype_table.at[1,"EM: A"])
            else:
                HLA_A_2 = "None"

            #B alleles
            HLA_B_1 = re.sub(pattern, "", hla_genotype_table.at[0,"EM: B"])
            if idx > 0:
                if hla_genotype_table.at[1, "EM: B"] == "None":
                    HLA_B_2 = HLA_B_1
                else:
                    HLA_B_2 = re.sub(pattern, "", hla_genotype_table.at[1, "EM: B"])
            else:
                HLA_B_2 = "None"

            #C alleles
            HLA_C_1 = re.sub(pattern, "", hla_genotype_table.at[0,"EM: C"])
            if idx > 0:
                if hla_genotype_table.at[1, "EM: C"] == "None":
                    HLA_C_2 = re.sub(pattern, "", hla_genotype_table.at[0, "EM: C"])
                else:
                    HLA_C_2 = re.sub(pattern, "", hla_genotype_table.at[1, "EM: C"])
            else:
                HLA_C_2 = "None"

        #get the KIR alleles
        #if kir_genotype_table.loc["KIR2DL1", "qscore_1"] > 0:
        KIR2DL1_1 = kir_genotype_table.loc["KIR2DL1","allele_1"]
        KIR2DL1_2 = kir_genotype_table.loc["KIR2DL1","allele_2"]
        KIR2DL2_1 = kir_genotype_table.loc["KIR2DL2","allele_1"]
        KIR2DL2_2 = kir_genotype_table.loc["KIR2DL2","allele_2"]
        KIR2DL3_1 = kir_genotype_table.loc["KIR2DL3","allele_1"]
        KIR2DL3_2 = kir_genotype_table.loc["KIR2DL3","allele_2"]
        KIR3DL1_1 = kir_genotype_table.loc["KIR3DL1","allele_1"]
        KIR3DL1_2 = kir_genotype_table.loc["KIR3DL1","allele_2"]

        new_row = {"sequencing_id":sequencing_id, "HLA-A_1":HLA_A_1, "HLA-A_2":HLA_A_2, "HLA-B_1":HLA_B_1, "HLA-B_2":HLA_B_2, "HLA-C_1":HLA_C_1, "HLA-C_2":HLA_C_2, "KIR2DL1_1":KIR2DL1_1, "KIR2DL1_2":KIR2DL1_2, "KIR2DL2_1":KIR2DL2_1, "KIR2DL2_2":KIR2DL2_2, "KIR2DL3_1":KIR2DL3_1, "KIR2DL3_2":KIR2DL3_2, "KIR3DL1_1":KIR3DL1_1, "KIR3DL1_2":KIR3DL1_2}    

    return new_row

#create empty dataframe to store the patient alleles, each row is a patient
patient_alleles = pandas.DataFrame(columns = ["sequencing_id","HLA-A_1","HLA-A_2","HLA-B_1", "HLA-B_2","HLA-C_1","HLA-C_2","KIR2DL1_1","KIR2DL1_2","KIR2DL2_1","KIR2DL2_2","KIR2DL3_1","KIR2DL3_2","KIR3DL1_1","KIR3DL1_2"])

# #iterate through the rows of the dataframe, each row has a pair of files, run the pair through the function
for index, (patient_hla_file, patient_kir_file) in files.iterrows():
    #patient row has all of a patient's relevant alleles
    new_row = extract_patient_alleles(patient_hla_file, patient_kir_file)
    #add the row to the dataframe
    patient_alleles = pandas.concat([pandas.DataFrame([new_row]), patient_alleles], ignore_index = True)

#patient_alleles.to_csv("~/patient_alleles.csv")

#extract the sequences for the alleles present in the patient dataset
def get_allele_sequences(allele, sequence_file):
    for index, row in patient_alleles.iterrows():
        allele_name = row[allele].removeprefix("HLA-")
        for sequence in sequence_file:
            sequence_name = sequence[0]
            #if the patient allele matches an allele name with sequence
            if allele_name in sequence_name: 
                #add the sequence and allele to a new list
                patient_alleles_and_sequence.append([allele_name, sequence[1]]) 

patient_alleles_and_sequence = []

for column in patient_alleles:
    if "HLA-A" in column:
        get_allele_sequences(column, HLA_A_sequences)
    if "HLA-B" in column:
        get_allele_sequences(column, HLA_B_sequences)
    if "HLA-C" in column:
        get_allele_sequences(column, HLA_C_sequences)


#create empty columns in the patient dataframe to store the count and score
patient_alleles["motif"] = ""
patient_alleles["functional_ikir_count"] = 0
patient_alleles["score"] = 0

alleles_and_motifs = []

for allele_sequence in patient_alleles_and_sequence:
    allele = allele_sequence[0]
    sequence = allele_sequence[1]
    if len(sequence) > 103:
        if "SLRN" in sequence:
            if "C" in allele:
                alleles_and_motifs.append([allele, "C1"])
            if "B*46" in allele or "B*73" in allele:
                alleles_and_motifs.append([allele, "C1"])
        if "NLRK" in sequence and "C" in allele:
            alleles_and_motifs.append([allele, "C2"])
        if sequence[103] == "T" or sequence[103] == "I":
            if "B" in allele:
                alleles_and_motifs.append([allele, "Bw4"])
            if "A*23" in allele or "A*24" in allele or "A*32" in allele:
                alleles_and_motifs.append([allele, "Bw4"])
        elif sequence[79] == "T" or sequence[79] == "I":
            if "B" in allele:
                alleles_and_motifs.append([allele, "Bw4"])
            if "A*23" in allele or "A*24" in allele or "A*32" in allele:
                alleles_and_motifs.append([allele, "Bw4"])

allele_motif_dict = {allele: motif for allele, motif in alleles_and_motifs}


# Iterate over the rows of the DataFrame
for index, row in patient_alleles.iterrows():
    if "KIR2DL1" in row["KIR2DL1_1"]:  #check if 'KIR2DL1' present
        #iterate over the alleles in the row
        for col in ["HLA-C_1", "HLA-C_2"]:
            allele = row[col]
            # If the allele is in the dictionary and the corresponding motif is 'C2', set 'motif' to 'C2'
            if allele in allele_motif_dict and allele_motif_dict[allele] == 'C2':
                patient_alleles.at[index, 'motif'] += 'C2 '
                patient_alleles.at[index, 'functional_ikir_count'] += 1
                patient_alleles.at[index, 'score'] += 1
                break  #exit the loop once the 'C2' motif is found
    if "KIR2DL2" in row["KIR2DL2_1"]:  #check if 'KIR2DL1' prsent
        #iterate over the alleles in the row
        for column in ["HLA-A_1", "HLA-A_2", "HLA-B_1", "HLA-B_2", "HLA-C_1", "HLA-C_2", "motif"]:
            allele = row[column]
            #if the allele is in the dictionary and the corresponding motif is 'C2', add a placeholder 
            #because C2/KIR2DL2 is the weak pairing, we'll add the score and count later if C1 is not present
            if allele in allele_motif_dict and allele_motif_dict[allele] == 'C2':
                patient_alleles.at[index, 'motif'] += 'C2 ' #we are going to add this as a placeholder so that it still searches for allele with C1, so no adding score/count yet
            #if the allele is in the dictionary and the corresponding motif is 'C1', then add the score and pair 
            elif allele in allele_motif_dict and allele_motif_dict[allele] == 'C1': #if it finds the C1, then add the pair and score
                patient_alleles.at[index, 'motif'] += 'C1 '
                patient_alleles.at[index, 'functional_ikir_count'] += 1 
                patient_alleles.at[index, 'score'] += 1 #this is the strong pairing
                break
    if "KIR2DL3" in row["KIR2DL3_1"]:  #check if 'KIR2DL3' present
            for col in ["HLA-A_1", "HLA-A_2", "HLA-B_1", "HLA-B_2", "HLA-C_1", "HLA-C_2"]:
                allele = row[col]
                #if the allele is in the dictionary and the corresponding motif is 'C2', set 'motif' to 'C2'
                if allele in allele_motif_dict and allele_motif_dict[allele] == 'C1':
                    patient_alleles.at[index, 'motif'] += 'C1 '
                    patient_alleles.at[index, 'functional_ikir_count'] += 1 
                    patient_alleles.at[index, 'score'] += 0.75
                    break
    if "KIR3DL1" in row["KIR3DL1_1"]:  #check if 'KIR2DL1' prsent
            #iterate over the alleles in the row
            for col in ["HLA-A_1", "HLA-A_2", "HLA-B_1", "HLA-B_2", "HLA-C_1", "HLA-C_2"]:
                allele = row[col]
                # If the allele is in the dictionary and the corresponding motif is 'C2', set 'motif' to 'C2'
                if allele in allele_motif_dict and allele_motif_dict[allele] == 'Bw4':
                    patient_alleles.at[index, 'motif'] += 'Bw4 '
                    patient_alleles.at[index, 'functional_ikir_count'] += 1 
                    patient_alleles.at[index, 'score'] += 1
                    break

for index, row in patient_alleles.iterrows():
    if "KIR2DL2" in row["KIR2DL2_1"]:  #check if 'KIR2DL2' prsent
        for column in ["motif"]: #now we are adding the score/count based on whether an allele with a C1 motif is present
            motif = row[column]
            if "C2" in motif and "C1" not in motif: #if the placeholder is there and C1 is not, THEN add the score for the weaker pairing
                patient_alleles.at[index, 'functional_ikir_count'] += 1
                patient_alleles.at[index, 'score'] += 0.5

print(patient_alleles)

#put them in a CSV to look at for now
#patient_alleles.to_csv(r"~/patient_alleles.csv")