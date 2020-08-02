# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | read1 |
| 1294_S1_L008_R2_001.fastq.gz | index1 |
| 1294_S1_L008_R3_001.fastq.gz | index2 |
| 1294_S1_L008_R4_001.fastq.gz | read2 |

2. Per-base NT distribution
    1. A per base distribution of quality scores for read1, read2, index1, and index2 are show below.
    
        read1
        ![](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-Chris-Anthony/master/Assignment-the-first/Mean_Qscore_R1.png)
        
        index1
        ![](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-Chris-Anthony/master/Assignment-the-first/Mean_Qscore_R2.png)
        
        index2
        ![](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-Chris-Anthony/master/Assignment-the-first/Mean_Qscore_R3.png)
        
        read2
        ![](https://raw.githubusercontent.com/2020-bgmp/demultiplexing-Chris-Anthony/master/Assignment-the-first/Mean_Qscore_R4.png)
    2. A good quality score cutoff for the index reads and biological read pairs to utilize for sample identification and downstream analysis for both would be Q30. A Q score of 30 is equivalent to the probability of an incorrect base call of 1 in 1000 times. A lower qscore (Q20) will have an incorrect base call probability of 1 in 100. For biological reads greater than 100 bp, there will be one inccorect baseEach base for all index and biological read pairs have a mean quality score above 30. However, in the early cycles, the mean quality score are close to this cutoff. In order to not remove these reads, I will take the average quality of the index reads.
    3. The number of indexes that have an undetermined (N) base calls are listed below along with the unix command.
    
    ```
    ls -1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R[2-3]_001.fastq.gz | while read file; do echo $file; zcat $file | sed -n "2~4p" | grep "N" | wc -l; done
    ```
    /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
    3976613
    /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
    3328051
    
## Part 2
1. Define the problem
2. Describe output
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [4 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
