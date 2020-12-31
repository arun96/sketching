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

# Java Implementation
There are three distinct approaches used to generate the screen - a MinHash-based approach, a Minimizer-based approach, and a Uniform sampling approach. For the first two, we have the option of calculating the sketch/window size, but also the option to specify the sketch/window size that will be used.

Now, there are multiple options for what to do with the generated screen.

-  To just generate a screen, and not classify any reads, please also see the "Screen/Sketch Generation Only" section of the README.
- To run the kinds of experiments I ran, select a set of genomes, simulate reads as specified above, and follow the instructions below to test the classification of those reads.
- To classify a brand new set of reads against this screen, see the "Classifying novel reads" section towards the end of the README.

As always, make sure to compile: `javac -cp jars/\* screen_java/*.java`.

For details on how to control the number of threads to use, look for the "Threads" section at the end of this README.

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
java -cp screen_java:jars/\* Main <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> <Hash Function>
```

### MinHash with specified sketch size

To use a MinHash-based approach to screen a generated set of reads against the set of genomes they were generated from, but with a specified sketch size for all genomes, use the following syntax:
```
java -cp screen_java:jars/\* Main <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> f <Fixed Sketch Size> <Hash Function>
```

### Minimizer with calculated window size

To use a Minimizer-based approach to screen a generated set of reads against the set of genomes they were generated from, use the following syntax:
```
java -cp screen_java:jars/\* Main <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> m <Hash Function>
```

### Minimizer with specified window size

To use a Minimizer-based approach to screen a generated set of reads against the set of genomes they were generated from, but with a specified window size, use the following syntax:
```
java -cp screen_java:jars/\* Main <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> m <Specified Window Size> <Hash Function>
```

### Uniform Sampling

To use a Uniform Sampling-based approach to screen a generated set of reads against the set of genomes they were generated from, use the following syntax:
```
java -cp screen_java:jars/\* Main <Genomes Directory/> <File to Save Output to> <Reads Directory/> <Read Length> <Read Error Rate> <Number of Target Matches per Read> u <Hash Function>
```

# ADDITIONAL OPTIONS

## Threads

By default, this program using 4 threads for read classification. I am working on making this a run-time parameter, but until then, you can change the number of threads by editing Line 85 in `Settings.java`. I will update this README when this changes!

## Loading fixed number of reads

There is an option to load reads in fixed-size chunks, instead of a whole file at a time. The relevant lines are in `Settings.java`, lines 89-91 - set `IN_CHUNKS = true` to enable this option, then set `CHUNK` to be the number of reads to be loaded and processed at once, and finally `CHUNK_UPDATES = true` if you want an update on the classification rate to be printed after each chunk (if `false`, then the classification rate will only be printed after the entire read set is processed). By default, the reads will be loaded and classified in chunks of 2000.

## Screen/Sketch Generation Only

If you wish to only generate screens/sketches for selected genomes, but not compare reads against them, use the following command:
```
java -cp mashscreen_java:jars/\* Main S <Genomes Directory/> <Folder to save screens to/> <Expected Read Length> <Expected Read Error Rate> <Expected Number of Target Matches per Read> <m/u/BLANK> <Hash Function to use>
```

The screens/sketches will be saved in the specified directory in `.bin` files. I am currently working on adding functionality for loading saved screens, so that will be added very soon.

## Read Logging

To save results for each individual read, please modify the parameters `READ_LOGGING` and `READ_LOCATION` in `Setting.java`. If the latter is set to `true`, then for each read, a `.log` file will be saved with the true source of the read, the predicted source, and the number of matches in each case. By default, `READ_LOGGING` is set to `false` if the read sets are generated from the genomes used in the screen (and therefore a ground truth is known), but log files will be saved to `READ_LOCATION` if new read sets are being classified (with no ground truth).

## Classifying Novel Reads

If the number of read sets does not match the number of elements in the screen, then we treat these read sets as novel, and therefore cannot be sure which genomes they truly come from. The code will still classify every read in every readset, but instead of outputting results for the number of classified/misclassified reads, it will save the output of each read's classification (i.e. how many matches it has to each element of the screen) to the location specified in `READ_LOCATION` in `Setting.java`. These results can be analyzed using the provided Python helper scripts (outlined in the "Analyzing Results" section below).

## Using External JARs

As I am not using Maven, I have manually included the JARs this project will use. The foremost of these is Google's `Guava`, which gives me access to hash function implementations, optimized data structures, and other nifty features. As I add more to this folder, I will update the README to include descriptions of each of them!

# Analyzing Results

TODO - this is still new, and I will update this over the coming weeks.

In the `analysis` folder, you can find some scripts that are useful for analyzing the results of the screening process.
- `extract_results.py`: Use `python3 extract_results.py <Path to output file>` to generate a summary of the experiment. This will give you the total number of correctly and incorrectly classified reads, as well as a few other metrics. It will also generate a histogram with the classification accuracy of reads from each organism.
- `aggregate_read_logs.py`: Use `python3 aggregate_read_logs.py <Path to folder with read logs> <numbers of members in the community` to get a matrix showing how many reads were correctly or incorrectly classified to each member of the community. This is only available if you generated individual log files for each read, by setting the `READ_LOGGING` option in `Settings.java` to `true`. The matrix is interpreted as follows - the value at` Matrix[x][y]` is the number of reads from organism `x` that were classified as being from organism `y`.
- `aggregate_classification_logs.py`: Use `python3 aggregate_classification_logs.py <Path to folder with read logs>` to generate a series of histograms showing the breakdown of how many reads were classified to each genome in the screen for each readset that was classified. The folder containing these log files will be the location stored in the `READ_LOCATION` parameter in `Settings.java`.
