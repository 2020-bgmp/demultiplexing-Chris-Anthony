#!/usr/bin/env python

print("Initializing...")
print()

def get_args():
    '''Argparse code to allow for in-line command interface.'''
    import argparse

    parser = argparse.ArgumentParser(description="Takes read1, read2, read3, and read4 FASTQ files. Returns analysis of quality scores by cycle.")
    parser.add_argument('-r','--read',  action='store', nargs='?', type=str, 
                    required=True, help='read file that ends in .fastq', dest="read")

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


def displayResults(mean_scores:list)->None:
    '''Takes a list of of mean qscores (list) per FASTQ file 
    and prints the results and distribution.'''

    # Print results
    print("# Base", "Mean Qscore", sep="\t")

    for i in range(len(mean_scores)):
        print(i, format(mean_scores[i],'.2f'), sep='\t')
    
    # Plot distribution
    import matplotlib.pyplot as plt

    plot1 = plt.figure(1)
    plt.plot(mean_scores)
    plt.title("Average Quality Score Per Base")
    plt.ylabel("Mean Quality Score")
    plt.xlabel("# Base")
    plt.savefig("Mean_Qscore")

    return None


def main(file1:str):
    '''Takes read1, read2, read3, and read4 FASTQ files. 
    Returns analysis of quality scores by cycle.'''

    l = seqLength(file1)
    list1, LN = populate_list(file1, l)
    mean_scores = statList(list1, LN, l)
    displayResults(mean_scores)


if __name__ =="__main__":
    args = get_args()
    read = args.read

    main(read)