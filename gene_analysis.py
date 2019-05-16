# gene_analysis.py
# Author: Eden Johnson
# Last updated: 5/15/19
# Purpose: Program references the Alzheimer's Disease gene coding sequence and
#          finds the % of each base in the gene, counts all the dinucleotides,
#          counts of the trinucleotides, the AT content of the gene and
#          presence of any repeats to provide analysis. 
# Program uses: Built-in Python functions open(), len(), print(), str(),
#               round() built-in Python methods, .readlines(), .rstrip(),
#               .write(), .read(), .upper(), .count(). Use of dictionarries,
#               conditional statement, for loops. Import of the re module

import re

# open the gene file and skip over the 1st line of comments
AD_file = open("AD_sequence.fasta", 'r')
AD_coding = open("AD_coding.txt", 'w')
file = AD_file.readlines()[1:]
# write the coding sequence to a new file
for line in file:
    line_data = line.rstrip("\n")
    AD_coding.write(line)

# close the isolated DNA file
AD_coding.close()

# reopen the file for reading
dna_file = open("AD_coding.txt", 'r')
dna = dna_file.read()

# refer a base list to count dinucleotides and trinucleotides
base_list = ['A','T','G','C']
# Initialize the dictionary all_counts, di_counts and tri_counts
di_counts = {}
tri_counts = {}
all_counts = {}

# get base percentage 
def get_base_perecent(dna,base):
    """ Function gets the percentage of bases in the strad of DNA. """
    length = len(dna)
    count = dna.upper().count(base)
    percentage = count/length*100
    return "%.3f" % percentage
for base in base_list:
    print("The percentage of " + base + " is: " + get_base_percent(dna,base)
          + "%")

# count the dinucleotides
for base_1 in base_list:
    for base_2 in base_list:
        dinucleotide = base_1 + base_2
        count = dna.count(dinucleotide)
        if count > 0:
            all_counts[dinucleotide] = count
        # assign value to di dictionary key
        di_counts[dinucleotide] = count
for di in di_counts:
    print("The count for " + di + " is: " + str(all_counts.get(di,0)))

# count the trinucleotides    
for base_1 in base_list:
    for base_2 in base_list:
        for base_3 in base_list:
            # create a trinucleotide
            trinucleotide = base_1 + base_2 + base_3
            count = dna.count(trinucleotide)
            if count > 0:
                all_counts[trinucleotide] = count
            # assign value to tri dictionary key
            tri_counts[trinucleotide] = count
for tri in tri_counts:
    print("The count for " + tri + " is: " + str(all_counts.get(tri,0)))

# get the AT content of the gene
def get_at_content(dna):
    """Function computes the proportion of A's and T's in dna. """
    # body of the function
    length = len(dna)
    a_count = dna.upper().count("A")
    t_count = dna.upper().count("T")
    at_content = (a_count + t_count)/length*100.0
    # at_content is what this function outputs
    return "%.3f" % at_content 
print("The AT content of the gene is: " + str(get_at_content(dna)) + "%")

# find any repeats in the DNA
a_repeat_list = re.finditer(r"A{6,}",dna)
t_repeat_list = re.finditer(r"T{6,}",dna)
c_repeat_list = re.finditer(r"C{6,}",dna)
g_repeat_list = re.finditer(r"G{6,}",dna)

# check for each base for repeats
for repeat in a_repeat_list:
    print("Repeat of " + repeat.group() + " found at position: " +
          str(repeat.start()))
for repeat in t_repeat_list:
    print("Repeat of " + repeat.group() + " found at position: " +
          str(repeat.start()))
for repeat in c_repeat_list:
    print("Repeat of " + repeat.group() + " found at position: " +
          str(repeat.start()))
for repeat in g_repeat_list:
    print("Repeat of " + repeat.group() + " found at position: " +
          str(repeat.start()))

# close the files
AD_file.close()
dna_file.close()
