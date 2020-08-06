#!/usr/bin/env python

print()
print("Initializing...")

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

def createBarcodePairs(indexes:str)->(dict, list):
    '''Take in the indexes.txt, returns a list of the indexes and a list of samples'''
    import re
    import itertools

    DNA = re.compile("[ATCG]{2,}")

    temp = []
    index_list = []
    sample_dict = {}

    with open(indexes, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            temp = line.split("\t")
            if DNA.search(line):
                sample_dict[temp[4]] = temp[0]
                index_list.append(temp[4])
    
    index_list = list(itertools.product(index_list, repeat=2))

    for i in range(len(index_list)):
        index_list[i] = index_list[i][0] + "-" + index_list[i][1]

    return sample_dict, index_list

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
    '''takes barcode records and attempts to error correct'''

    # Set froward and reverse equal to index 1 
    # and reverse complement of index 2, respectively
    # remove new line character
    forward = record[1].strip()
    reverse = reverseComp(record[2].strip())

    for i in range(0, len(forward)):
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
    '''Takes the record from a list of FASTQ files (2 index reads and 2 biological reads).
        Demultiplexes record whether barcodes match/mismatch and quality is above/below the threshold.'''
    
    # If there are N in the indexes, attempt to error correct
    if "N" in record[1][1] or "N" in record[1][2]:
        record[1] = errorCorrect(record[1])

    # Call function to add barcode pair to header lines
    record, barcode_pair = add2header(record)
    
    # Count number of occurances per barcode pair
    check =[]

    if barcode_pair in dict_count.keys():
        dict_count[barcode_pair][0] += 1
        check = barcode_pair.split("-")
        if check[0] != check[1]: # count number of mismatches
            barcode_pair = "mismatch"
            dict_count["mismatch"][0] += 1
    else: # for barcode pairs not in dictionary (e.g. contains N), count as unknown
        barcode_pair = "lowQ_unk"
        dict_count["lowQ_unk"][0] += 1
    
    # Write to output files
    if meanQual(record[3][1].strip()) < qscore or meanQual(record[3][2].strip()) < qscore or barcode_pair == "lowQ_unk":
        barcode_pair = "lowQ_unk" # low quality or unknown
        dict_count = exportRecord(record, barcode_pair, dict_count)
    else: 
        if barcode_pair == "mismatch": # mismatches
            dict_count = exportRecord(record, barcode_pair, dict_count)
        else: # matches
            dict_count = exportRecord(record, barcode_pair, dict_count)

    return dict_count

def printAnalysis(dict_count:dict, sample_dict:dict)->None:
    '''Print a demultiplex report as a textfile using a sample list and dictionary of barcode pair occurances.'''
    # Count the number of matched indexes
    match = 0
    check = []
        
    for key in dict_count.keys():
        if key not in ["mismatch", "lowQ_unk"]:
            check = key.split("-")
            if check[0] == check[1]:
                match += dict_count[key][0]

    total = match + dict_count["mismatch"][0] + dict_count["lowQ_unk"][0]

    with open("demultiplex_results.txt", "w") as OUT:
        OUT.writelines(["The number of read-pairs with matched indexes are (% total reads):\t", str(match), "\t(", str(format(match/total, '.2%')), ")\n"])
        OUT.writelines(["The number of read-pairs with mismatched indexes are (% total reads):\t", str(dict_count["mismatch"][0]), "\t(", str(format(dict_count["mismatch"][0]/total, '.2%')), ")\n"])
        OUT.writelines(["The number of read-pairs with unknown indexes after error correction are (% total reads):\t", str(dict_count["lowQ_unk"][0]), "\t(", str(format(dict_count["lowQ_unk"][0]/total, '.2%')), ")\n"])
        OUT.write("NOTE: Low quality reads were not separated from these counts.\n\n")

        OUT.writelines(["Sample\t", "Barcode Pair\t", "Number of Occurances\t", "Percentage of Total\t","\n"])
        for key in sample_dict.keys():
            barcode_pair = key + "-" + key
            OUT.writelines([str(sample_dict[key]), "\t", str(barcode_pair), "\t", str(dict_count[barcode_pair][0]), "\t", str(format(dict_count[barcode_pair][0]/total, '.2%')), "\n"])

        OUT.writelines(["\nBarcode Pair\t", "Number of Occurances\t", "Percentage of Total\t","\n"])
        for key in dict_count.keys():
            if key != "mismatch":
                if  key == "lowQ_unk":
                    OUT.writelines(["unknown barcodes\t", str(dict_count[key][0]), "\t", str(format(dict_count[key][0]/total, '.2%')), "\n"])
                else:
                    OUT.writelines([key, "\t", str(dict_count[key][0]), "\t", str(format(dict_count[key][0]/total, '.2%')), "\n"])

    return None

def main(read_list:list, indexes:str, qscore:int):
    '''Takes in four read FASTQ files as a list, iterates record by record, 
    and separates records into output files based on barcode pair/quality'''
    import gzip
    import itertools

    # initialize dictionary
    dict_count = {}

    sample_dict, index_list = createBarcodePairs(indexes)

    for index in index_list:
        dict_count[index] = [0, "output_R1_" + index + ".fastq", "output_R4_" + index + ".fastq"]
    dict_count["lowQ_unk"] = [0, "output_R1_" + "lowQ_unk" + ".fastq", "output_R4_" + "lowQ_unk" + ".fastq"]
    dict_count["mismatch"] = [0, "output_R1_" + "mismatch" + ".fastq", "output_R4_" + "mismatch" + ".fastq"]

    # open FASTQ files simultaneously as store each record as a list
    record = []

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
    printAnalysis(dict_count, sample_dict)

    # Close all files
    for fh in files:
        fh.close()
    
    for key in dict_count:
        if type(dict_count[key][1]) != str:
            dict_count[key][1].close()
            dict_count[key][2].close()

    print("Complete.")

if __name__ =="__main__":
    args = get_args()
    read_list = [args.read1, args.read2, args.read3, args.read4]
    indexes = args.indexes
    qscore = args.qscore

    main(read_list, indexes, qscore)