#!/usr/bin/env python

print("Initializing...")
print()

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="""Takes read1, read2, read3, and read4 FASTQ files and indexes textfile. 
                    Demultiplexes textfiles and returns report of the match, mismatched and low quality/unknown barcodes""")
    parser.add_argument('-r1','--read1',  action='store', nargs='?', type=str, 
                    required=True, help='read1 file that ends in .fastq', dest="read1")
    parser.add_argument('-r2','--read2',  action='store', nargs='?', type=str, 
                    required=True, help='read2 file that ends in .fastq', dest="read2")
    parser.add_argument('-r3','--read3',  action='store', nargs='?', type=str, 
                    required=True, help='read3 file that ends in .fastq', dest="read3")
    parser.add_argument('-r4','--read4',  action='store', nargs='?', type=str, 
                    required=True, help='read4 file that ends in .fastq', dest="read4")
    parser.add_argument('-i','--index',  action='store', nargs='?', type=str, 
                    required=True, help='index textfile', dest="indexes")
    parser.add_argument('-q','--qscore',  action='store', nargs='?', type=int, 
                    default=30, help='Qscore threshold, default is Q30', dest="qscore")

    return parser.parse_args()

def reverseComp(seq:str)->str:
    '''Take in a sequence string. Returns its reverse complement.'''
    
    # Reverse the sequence
    reverse = seq[::-1]
    
    # Watson-Crick base pair as a dictionary
    convert_dict = {
        "A":"T",
        "T":"A",
        "C":"G",
        "G":"C",
        "N":"N"
        }

    # Watson-Crick base pair conversion 
    for i in range(len(reverse)):
        if i == 0:
            reverseComp = convert_dict[reverse[i]]
        else:
            reverseComp += convert_dict[reverse[i]]

    return reverseComp

def createBarcodePairs(indexes:str)->list: #Defunct, remove?
    '''Take in the indexes.txt, returns a list of the indexes'''
    import re

    DNA = re.compile("[ATCG]{2,}")

    index_list = []

    with open(indexes, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            index = DNA.findall(line)
            if index != []:
                index_list.append(index[0])
    
    return index_list

def convert_phred(letter:str)->int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33

def meanQual(line:str)->float:
    '''Takes a quality line and returns the average qvalue'''
    total = 0.0

    for ch in line:
        total += convert_phred(ch)
    
    return total/len(line)

def errorCorrect(record:list)->list:
    '''takes barcode records and attempts to eror correct'''

    # Set froward and reverse equal to index 1 
    # and reverse complement of index 2, respectively
    # remove new line character
    forward = record[1].strip()
    reverse = reverseComp(record[2].strip())

    for i in range(len(forward)):
        if i == 0: # Initializes new string
            if forward[i] != reverse[i]: 
                if forward[i] == "N": # If forward string starts with N, replace with reverse
                    seq = reverse[i]
                if reverse[i] == "N": # and vice versa
                    seq = forward[i]
            else:
                seq = forward[i] 
        else: # Concatenate to string while error correcting
            if forward[i] != reverse[i]:
                if forward[i] == "N":
                    seq += reverse[i]
                if reverse[i] == "N":
                    seq += forward[i]
            else:
                seq += forward[i]
    
    # Replace the original records with the error corrected sequences
    record[1] = seq + "\n"
    record[2] = reverseComp(seq) + "\n"

    return record

def add2header(record:list)->(list, str):
    '''Takes a record stored in a list. Adds barcode pair to header line.
        Returns barcode pair.'''

    # Create a barcode pair where indexes are separated by a dash
    barcode_pair = record[1][1].strip() + "-" + reverseComp(record[1][2].strip())

    # Add barcode pairs to header lines
    record[0][0] = record[0][0].strip() + ":" + barcode_pair + "\n"
    record[0][3] = record[0][3].strip() + ":" + barcode_pair + "\n"

    return record, barcode_pair

def exportRecord(record:list, barcode_pair:str, dict_count:dict)->dict:
    '''Export record to FASTQ file determined by barcode pair.'''

    # If the file is a string, open the file name and replace with its file object
    # All files are closed at the end of the main() function
    if type(dict_count[barcode_pair][1]) == str:
        dict_count[barcode_pair][1] = open(dict_count[barcode_pair][1], "w")
        dict_count[barcode_pair][2] = open(dict_count[barcode_pair][2], "w")
    
    # Now, write the writes to the output file
    for i in range(len(record)):
        dict_count[barcode_pair][1].writelines(record[i][0])
        dict_count[barcode_pair][2].writelines(record[i][3])

    return dict_count

def analyzeRecord(record:list, qscore:int, dict_count:dict)->dict:
    '''Takes the record from a set of FASTQ files (2 index reads and 2 biological reads).
        Demultiplexes record whether barcodes match/mismatch and quality is above/below the threshold.'''
    
    # If there are N in the indexes, attempt to error correct
    if "N" in record[1][1] or "N" in record[1][2]:
        record[1] = errorCorrect(record[1])

    # Call function to add barcode pair to header lines
    record, barcode_pair = add2header(record)
    list1 =[]

    # Update dictionary containing barcode pairs, values are output files
    # for low quality of unknown indexes
    dict_count["lowQ_unk"] = [0, "output_R1_" + "lowQ_unk" + ".fastq", "output_R4_" + "lowQ_unk" + ".fastq"]
    if "N" in barcode_pair:
        barocde_pair = "lowQ_unk"
        dict_count["lowQ_unk"][0] += 1

    # for matched indexes
    list1 = barcode_pair.split("-")
    if list1[0] == list1[1]:
        if barcode_pair not in dict_count:
            dict_count[barcode_pair] = [0, "output_R1_" + barcode_pair + ".fastq", "output_R4_" + barcode_pair + ".fastq"]
            dict_count[barcode_pair][0] += 1 #retaining barcode_pair as key as it is used to specify outfile
        else:
            dict_count[barcode_pair][0] += 1
    # for mismatched indexes
    else:
        barcode_pair = "mismatch"
        if barcode_pair not in dict_count:
            dict_count["mismatch"] = [0, "output_R1_" + "mismatch" + ".fastq", "output_R4_" + "mismatch" + ".fastq"]
            dict_count["mismatch"][0] += 1
        else:
            dict_count["mismatch"][0] += 1
    
    # Write to output files
    # Check average quality of barcodes
    if meanQual(record[3][1].strip()) < qscore or meanQual(record[3][2].strip()) < qscore or barcode_pair == "lowQ_unk":
        barcode_pair = "lowQ_unk"
        dict_count = exportRecord(record, barcode_pair, dict_count)
    else: 
        if barcode_pair == "mismatch":
            dict_count = exportRecord(record, barcode_pair, dict_count)
        else:
            dict_count = exportRecord(record, barcode_pair, dict_count)

    return dict_count

def main(read_list:list, indexes:str, qscore:int):
    '''Takes in four read FASTQ files, iterates record by record, 
    and separates records into output files based on barcode pair/quality'''
    # index_list = createBarcodePairs(indexes) # defunct - possibly used to double check keys?

    import gzip
    import itertools

    record = []
    dict_count = {}

    # open FASTQ files simultaneously
    files = [gzip.open(fh,"rt") for fh in read_list]
    for line in itertools.zip_longest(*files):
        line = list(line)
        if line[0][0] != "@":
            record.append(line) # append record to list
        else:
            if record != []: # if record is not empty, call analyze Function
                dict_count = analyzeRecord(record, qscore, dict_count)
            record = []
            record.append(line) # append header line
    dict_count = analyzeRecord(record, qscore, dict_count)

    # Print number of matched, mismatched and unknown index pairs
    match = 0

    # Count the number of matched indexes
    for key in dict_count.keys():
        if key not in ["mismatch", "lowQ_unk"]:
            match += dict_count[key][0]

    print("The number of read-pairs with matched indexes are:\t", match)
    print("The number of read-pairs with mismatched indexes are:\t", dict_count["mismatch"][0])
    print("The number of read-pairs with unknown indexes after error correction are:\t", dict_count["lowQ_unk"][0])

    # Close all files
    for fh in files:
        fh.close()
    
    for key in dict_count:
        if type(dict_count[key][1]) != str:
            dict_count[key][1].close()
            dict_count[key][2].close()

if __name__ =="__main__":
    args = get_args()
    read_list = [args.read1, args.read2, args.read3, args.read4]
    indexes = args.indexes
    qscore = args.qscore

    main(read_list, indexes, qscore)