#!/usr/bin/env python

print("Initializing...")
print()

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="Takes read1, read2, read3, and read4 FASTQ files. Returns analysis of quality scores by cycle.")
    parser.add_argument('-r1','--read1',  action='store', nargs='?', type=str, 
                    required=True, help='read1 file that ends in .fastq', dest="read1")
    parser.add_argument('-r2','--read2',  action='store', nargs='?', type=str, 
                    required=True, help='read2 file that ends in .fastq', dest="read2")
    parser.add_argument('-r3','--read3',  action='store', nargs='?', type=str, 
                    required=True, help='read3 file that ends in .fastq', dest="read3")
    parser.add_argument('-r4','--read4',  action='store', nargs='?', type=str, 
                    required=True, help='read4 file that ends in .fastq', dest="read4")

    return parser.parse_args()


def init_list(my_list:list, l:int, value=0.0)->list:
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with l number of values of 0.0.'''
    
    for x in range(l):
        my_list.append(value)
        
    return my_list


def convert_phred(letter:str)->int:
    '''Converts a single character into a phred score'''
    return ord(letter)-33


def seqLength(file1:str)->int:
    '''This function takes a FASTQ file 
    and determines the length of the sequence'''
    import gzip

    fh = gzip.open(file1, "rt")

    for i, line in enumerate(fh):
        if i == 0:
            continue
        elif i == 1:
            line = line.rstrip("\n")
            l = len(line)
        else:
            break #I think this will prevent reading the whole file

    fh.close()

    return l


def populate_list(file1:str, l:int)->(float,int):
    '''This function takes a FASTQ file and returns the 
    summation of the phred values per base in a list 
    and the number of lines as an integer'''
    import gzip

    all_qscores = []

    for x in range(l):
        all_qscores.append([])

    LN = 0

    with gzip.open(file1,"rt") as fh:
        for line in fh:
            LN += 1
            line = line.rstrip("\n")
            if LN%4 == 0:
                for i in range(len(line)):
                    all_qscores[i].append(convert_phred(line[i]))
    
    return all_qscores, LN


def statList(all_qscores:list, LN:int, l:int)->None:
    '''Takes a list of lists of phred scores. Calculates the mean, variance,
    standard deviation, and median. Prints the values.'''
    # Calculate Mean
    mean_scores = []
    mean_scores = init_list(mean_scores, l)
    
    for i in range(len(all_qscores)):
        for j in range(len(all_qscores[i])):
            mean_scores[i] += all_qscores[i][j]
        mean_scores[i] = 4*mean_scores[i]/LN
    
    return mean_scores


def displayResults(mean_list:list)->None:
    '''Takes a list of of mean qscores (list) per FASTQ file 
    and prints the results and distribution.'''

    # Print results
    print("Mean Quality Scores")
    print("# Base", "R1", "R2", "R3", "R4", sep="\t")

    for i in range(len(mean_list[0])):
        print(i, end="\t")

        for j in range(len(mean_list)):
            try:
                print(format(mean_list[j][i],'.2f'), end="\t")
            except IndexError:
                print(end='\t')
        
        print()
    
    # Plot distribution
    import matplotlib.pyplot as plt

    for i in range(len(mean_list)):
        plot1 = plt.figure(i+1)
        plt.plot(mean_list[i])
        plt.title("Average Quality Score Per Base")
        plt.ylabel("Mean Quality Score")
        plt.xlabel("# Base")
        plt.savefig("Mean_Qscore_R%i" %int(i+1))

    return None


def main(read_list:list):
    '''Takes read1, read2, read3, and read4 FASTQ files. 
    Returns analysis of quality scores by cycle.'''
    mean_list = []
    x, y = 1, 1

    print("Input File Name", "Paired End Read or Index", sep='\t')

    for file1 in read_list:
        l = seqLength(file1)
        list1, LN = populate_list(file1, l)
        mean_scores = statList(list1, LN, l)
        mean_list.append(mean_scores)

        # Determine which files contain the indexes and which contain the paired end reads
        print(file1, end='\t')
        if l > 16: # Illumina barcodes will never be greater than 16
            print("read", x, sep='')
            x += 1
        if l <= 16:
            print("index", y, sep='')
            y += 1

    print()
    displayResults(mean_list)


if __name__ =="__main__":
    args = get_args()
    read1 = args.read1
    read2 = args.read2
    read3 = args.read3
    read4 = args.read4
    read_list = [read1, read2, read3, read4]

    main(read_list)