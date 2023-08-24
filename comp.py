import glob 
import numpy as np 
import pandas as pd 
import csv 
import random 

def add_phylum(fasta_path, edited_path, pfam_ID):
	path = fasta_path
	path_edited = edited_path

	fasta_list = glob.glob(path + '*.fa')

	for item in fasta_list: 
		filename = item 

		with open(filename, 'r') as f: 
			phylum = filename.split('-')[-1].split('.')[0]
			name = pfam_ID + '_selected_sequences-'
			newfile = path_edited + name + phylum + '_edited' + filename[-3:]

			with open(newfile, 'w') as g: 
				for line in f: 
					if '>' in line: 
						g.write('>' + phylum + '_' + line[1:])
					else: 
						g.write(line)

	return print('completed')

#########################################

def sequence_counter(filename):
    count = 0
    with open(filename, 'r') as g: 
        for line in g: 
            if '>' in line: 
                count = count + 1
    return count 

########################################

def sequence_total(file_list):
    total = 0
    
    for file in file_list:   
        count = sequence_counter(file)
        total = total+count
    
    return total


#######################################

def pfam_merger(data_path, output_file):
    fasta_list = glob.glob(data_path + '*.fa')
    nf_name = output_file

    newfile = open(nf_name, 'w+')

    for item in fasta_list: 
        filename = item 
        
        with open(filename, 'r') as f: 
            with open(nf_name, 'a') as g: 
                for line in f: 
                    g.write(line)
    
    merged_counted = sequence_counter(nf_name)
    total_counted = sequence_total(fasta_list)
    
    if merged_counted == total_counted:
        return print('No data lost! Woohoo!')
    else: 
        return print('Something went wrong, RIP')

##############################################

def ncbi_merger(ncbi_path, pfam_path):
    ncbi_file = [ncbi_path]
    pfam_file = glob.glob(pfam_path+'*.fa')

    all_files = ncbi_file + pfam_file

    nf_name = 'AllFiles_merged.fa'

    newfile = open(nf_name, 'w+')

    for item in all_files: 
        filename = item 

        with open(filename, 'r') as f: 
            with open(nf_name, 'a') as g: 
                for line in f: 
                    g.write(line)
    
    ncbi_count = sequence_total(ncbi_file)
    pfam_count = sequence_total(pfam_file)
    total_count = sequence_counter(nf_name)
    
    if (ncbi_count + pfam_count) == total_count: 
        return print('No data lost! Woohoo!')
    else: 
        return print('Something went wrong, rip :(')

#################################################

def randomizer(cbiR_f, pfam_f, total_seq):     
    
    cbiR = pd.read_csv(cbiR_f)
    pfam = pd.read_csv(pfam_f)
    
    cbiR_list = cbiR['Accession'].values.tolist()
    pfam_list = pfam['Merged'].values.tolist()
    
    sum_list = cbiR_list + pfam_list
    
    random_list = random.sample(sum_list, k=total_seq)
    print(len(random_list))
    print(len(set(random_list)))
    
    return(random_list)

###################################################

def part_random(cbiR_f, pfam_f, cbiR_num, pfam_num):
    cbiR = pd.read_csv(cbiR_f)
    pfam = pd.read_csv(pfam_f)

    cbiR_list = cbiR['Accession'].values.tolist()
    pfam_list = pfam['Merged'].values.tolist()

    cbiR_rand = random.sample(cbiR_list, k=cbiR_num)
    pfam_rand = random.sample(pfam_list, k=pfam_num)

    final = cbiR_rand + pfam_rand

    return(final)
    
###################################################

def selector(input_file, output_file, input_list):
    with open(output_file, 'w+') as f: 
        with open(input_file, 'r') as g:
            header_set = set(input_list)
            i = 0
            lines = g.readlines()
            header_count = 0
            while i < len(lines):
                line = lines[i]
                if '>' in line and len(header_set) > 0:
                    header = next((ele for ele in header_set if ele in line), None)
                    if header is not None:
                        # write the header line if it's a header, otherwise just advance i
                        f.writelines(line)
                        # header_set.remove(ele)
                        header_count += 1
                        i += 1
                        while i < len(lines) and ('>' not in lines[i]):
                            f.writelines(lines[i])
                            i += 1
                    else:
                        i += 1
                else:
                    i += 1
            print(header_count)

####################################################

def x_reference(input_file, input_list):
    with open(input_file, 'r') as g: 
        for header in input_list:
            for line in g:
                if any(ele in line for ele in input_list):
                    continue
                else: 
                    print(header)


def y_reference(input_file, input_list):
    empty = []

    with open(input_file, 'r') as g: 
        for header in input_list: 
            count = 0
            for line in g: 
                if header in line: 
                    count += 1
                else: 
                    continue

            if count == 0:
                empty.append(header)
            else: 
                continue

    return empty


def missed(ref_list, cbiR_meta, pfam_meta):
    with open(input_file, 'r') as g: 
        cbiR = pd.read_csv(cbiR_meta)
        pfam = pd.read_csv(pfam_meta)

        cbiR_list = cbiR['Accession'].values.tolist()
        pfam_list = pfam['Pfam'].values.tolist()
        sum_list = cbiR_list + pfam_list

        diff = [x for x in sum_list if x not in ref_list]

    return print(diff)



####################################################

def pf_meta(dpath, output_file):
    
    fasta_list = glob.glob(dpath+'*.fa')
    
    with open(output_file, 'w', newline='') as csvfile: 
            header = ['Phylum', 'Pfam','Predicted Function', 'EMBL']
            writer = csv.DictWriter(csvfile, fieldnames=header, delimiter=',')
            writer.writeheader()
            
    for file in fasta_list: 
        with open(file, 'r') as f:
            with open(output_file, 'a', newline='') as csvfile: 
                writer = csv.DictWriter(csvfile, fieldnames=header, delimiter=',')

                for line in f: 
                    if '>' in line: 
                        line = line[1:]
                        bits = line.split()

                        pfam = bits[0].split('_')[-1]
                        protein = ' '.join(bits[1:]).split('{')[0]
                        phylum = file.split('-')[-1].split('_')[0] 
                        embl = line.split('{')[-1].split('}')[0]


                        data = {'Phylum': phylum, 'Pfam':pfam, 'Predicted Function': protein, 'EMBL':embl}
                        writer.writerow(data)
    
    check = pd.read_csv(output_file)
    
    if check['Phylum'].nunique() == len(fasta_list) and check['Pfam'].nunique() == sequence_total(fasta_list):
        return print('Everything looks ok!')
    else: 
        return print('Oops something is not ok :(')

#####################################################

def ncbi_meta(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

        count = 0 

        with open(output_file, 'a', newline='') as meta: 
            header = ['Accession', 'Protein', 'Species']
            writer = csv.DictWriter(meta, fieldnames=header, delimiter=',')
            writer.writeheader()

            for line in lines:
                if '>' in line: 
                    line = line[1:]
                    parts = line.split()

                    accession = parts[0]
                    species = line.split('[')[-1].split(']')[0]
                    protein = ' '.join(parts[1:]).split('[')[0] 

                    ncbi_dict = {'Accession':accession, 'Protein': protein, 'Species':species}
                    writer.writerow(ncbi_dict)

              
                    count = count + 1

            if count == sequence_counter(input_file):
                return print('No data lost! Yay!')
            else: 
                return print('Something went wrong :o')


####################################################





















