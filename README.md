# Overview
Code used to run the read classification experiments.

This is very much unfinished/unpolished, and I will continue to add to it and update the README to reflect the changes.

# Poster

The poster displaying results of this work can be found here: https://drive.google.com/file/d/1gFN_F26f3UZwSZawbVyWXuP3DnGB9_sH/view?usp=sharing

# Genomes
The MBARC-26, ZYMO, and MBARC + ZYMO genomes can be downloaded here: https://drive.google.com/drive/folders/1c-6B-G1-RGbIqzDhxzkV5smyv8XB32Xa?usp=sharing.

# Generating Reads
The read simulator is a modified version of Melanie Kirsche's read simulator (https://github.com/schatzlab/centroTools/tree/master/java). To generate reads for a set of genomes, run `generate_reads.sh` using the following syntax:
```
./generate_reads.sh <path to directory containing read simulator> <directory containing genomes> <directory where reads will be saved> <SNP rate> <Insertion rate> <Deletion rate> <Mean Read Length> <Coverage>
```

For example, `./generate_reads.sh readsim MBARC_ZYMO reads_MBARC_ZYMO 0.0034 0.0033 0.0033 10000 10` will simulate reads from all .fasta/.fna/.fa genome files in `MBARC_ZYMO`, and save a readset for each genome in `reads_MBARC_ZYMO`. The reads will have a mean length of 10K, a SNP rate of 0.34%, insertion and deletion rates of 0.33% each, and enough reads will be generated to have a coverage of 10x on each genome.

### Read Length Options (Not recommended)

By default, this script will generate reads with read lengths that are normally distributed around the input read length (with standard deviation equal to the square root of this mean length). However, you can generate reads with the following read length distributions:

- Exact [`XL`]: All reads will be the same length.
- Exponentially Distribution [`E`]: Read lengths will be drawn from an exponential distribution around the provided read length.
- Exponentially Distribution with Minimum Length [`EM`]: Read lengths will be drawn from an exponential distribution around the provided read length, with a minimum length of half the input length - if a read length is sampled that is lower than this minimum, that read's length will be the minimum value. This will result in many reads having read length exactly half the input length, so the next option is preferred for more realistic read lengths.
- Exponentially Distribution with Minimum Length v2 [`EL`]: Read lengths will be drawn from an exponential distribution around the provided read length, with a minimum length of half the input length. However, if a read length is below the minimum value, the length is sampled again. This creates more realistic read lengths.

By default, the normal distribution will be used. However, to use any of the other four options, simply add the appropriate string (`XL/E/EM/EL`) as an additional parameter when generating reads.

# JAVA
There are three distinct approaches used to generate the screen - a MinHash-based approach, a Minimizer-based approach, and a Uniform sampling approach. For the first two, we have the option of calculating the sketch/window size, but also the option to specify the sketch/window size that will be used.

First, make sure to compile: `javac screen_java/*.java`

#### Hash Function Options

At this point, the default hash function is Java's built in `hashCode()`. For the following section, if no hash function is specified, then this will be used. However, to specify a particular function, please use the following parameters:

- `h`: Java's `hashCode()`. Default option.
- `mmh3`: Google Guava's 32-bit `MurmurHash3` implementation.
- `mmh3_128`: Google Guava's 128-bit `MurmurHash3` implementation.

The selected hash function will be used to generate the screen and to classify the reads.

I will be adding more hash functions soon!

### MinHash with calculated sketch size

To use a MinHash-based approach to screen a generated set of reads against the set of genomes they were generated from, use the following syntax:
```
java -cp screen_java:jars/\* ReadScreener <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> <Hash Function>
```

### MinHash with specified sketch size

To use a MinHash-based approach to screen a generated set of reads against the set of genomes they were generated from, but with a specified sketch size for all genomes, use the following syntax:
```
java -cp screen_java:jars/\* ReadScreener <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> f <Fixed Sketch Size> <Hash Function>
```

### Minimizer with calculated window size

To use a Minimizer-based approach to screen a generated set of reads against the set of genomes they were generated from, use the following syntax:
```
java -cp screen_java:jars/\* ReadScreener <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> m <Hash Function>
```

### Minimizer with specified window size

To use a Minimizer-based approach to screen a generated set of reads against the set of genomes they were generated from, but with a specified window size, use the following syntax:
```
java -cp screen_java:jars/\* ReadScreener <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> m <Specified Window Size> <Hash Function>
```

### Uniform Sampling

To use a Uniform Sampling-based approach to screen a generated set of reads against the set of genomes they were generated from, use the following syntax:
```
java -cp screen_java:jars/\* ReadScreener <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> u <Hash Function>
```

# PYTHON
The Python implementation is not up to date, and should not be used.
