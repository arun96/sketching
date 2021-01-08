# Overview
Code used to run my sketching read classification experiments.

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
The Java implementation provides a way to generate a "screen" (sketched representation) of an input set of genomes, and either save the screen or classify input reads against this screen. For read classification, these may be reads drawn from the same set of genomes in the screen (either simulated using the approach above, or from another source), or a totally different set of reads.

There are three distinct approaches for screen-generation: a MinHash-based approach, a Minimizer-based approach, and a Uniform sampling approach. For each of these approaches, the size of the screen can either be calculated based on user-specified input parameters (the expected read length, the expected read error, the sizes of the genomes, and the number of shared samples each read should aim to have with its correct genome), or can be specified ahead of time by the user.

To see details of how to do any of the experiments outlined above, or how to adjust certain options/settings, please read the following sections of the README.

As always, make sure to compile: `javac -cp jars/\* screen_java/*.java`.

### Brief overview of the approach

For an input set of genomes, using the specified screen generation approach, a reduced representation of the genomes is generated. This may then be saved to a file.

If there are reads to be classified against the selected screen, the reads are loaded in (either one file at a time, or a specified number of reads at a time), and each read is compared against the reduced representations of the genome. The genome with which a given read shares the most hashed k-mers is the predicted source for the read.

If the reads and the genomes are "matched" (i.e. if there is one read set for each genome, and thus the ground truth is known), we are able to give a summary of the classification accuracy on each read set (we can also save a simple log file for each read, with the true source and the predicted source, and the number of matches with each). If the source of the reads is not known (i.e. they are "unmatched"), then a log file is generated for each read, showing the number of matches it had with each genome in the screen.

Finally, there are some Python scripts that can be used to aggregate and analyze individual read log files.

### Running the code and explanation of parameters

To run the code, please run:
```
java -cp screen_java:jars/\* Main <options>
```

This code uses the `commons.cli` argument parser, so the order of arguments does not matter. Below is a breakdown of the parameters/arguments that can be used.

Key Parameters:
- `-g/--genome <Directory>`: The directory containing the genomes that are to be sketched into a screen.
- `-r/--reads <Directory>`: The directory containing the reads that are to be classified against the screen. Only not necessary if we are in "screen-only" mode (`-so/--screen-only`).
- `-rl/--read-length <Integer>`: The expected read length of reads that will be classified against this screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-re/--read-error <Double>`: The expected error rate of reads that will be classified against this screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-tm/--target-matches <Integer>`: The target number of matches between a read and its correct source in the screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-s/--screen-type`: The method of screen-generation to be used for this screen. Choose between MinHash-based (default), [m]inimizer based, or [u]niform sampling.
- `-hf/--hash <Hash Type>`: The hash function to use throughout the process. Choose between Java's built in hashcode (default), [mmh3] (MurmurHash3) or [mmh3_128] (the 128-bit variant of MurmurHash3).

Matched Reads/Genomes vs Classification without a ground truth:
- `-um/--unmatched`: Use this parameter if the reads and the genomes do not correspond - by default, the code will assume that there is a matching read set for each genome (generated using the simulator above), and will proceed with classification assuming there is some ground truth. If this option is selected, detailed classification logs will be saved for each read, and no accuracy metrics will be printed. NOTE: Over time, this will become the default, and a parameter will be needed to indicate the genomes and reads are matched. In this mode, read logging is enabled by default.

Fixed Size Screens:
- `-f/--fixed <Integer>`: Use if the screen size should not be calculated, but instead the specified size should be used for all screens.

Read Loading:
- `-c/--chunks`: Use if a specified number of reads should loaded at a time (i.e. in "chunks"), instead of one file at a time.
- `-cs/--chunk-size <Integer>`: If `-c` is used, then `-cs` will specify the "chunk" size.
- `-cu/--chunk-updates`: If `-c` is used, then enabling `-cu` will print an update after each chunk is processed.

Read Logging:
- `-rlg/--read-logging`: Use if read-logging should be enabled. This is enabled by default in unmatched mode (`-um`).
- `rlc/--read-location <Directory>`: The directory that read logs should be saved to. They will be named with the naming convention `<Readset Number>_<Read_Number>.log` - for example, `3_23.log` is the log for the 24th read from the 3rd read set.

"Screen Only":
- `-so/--screen-only`: Use if you only want to generate the screen, but not classify any reads.
- `-sl/--screen-location <Directory>`: Specifies the location where the generated screens will be saved. Screens will be saved in `.bin` files with names matching the input genome files.
- `-cmbs/--combined-screens`: Save all screens as one large `.bin` file. By default, they are saved separately.

Experiment Parameters
- `-k/--kmer <Integer>`: The k-mer size to use. By default, it is 21.
- `-nt/--num-threads <Integer>`: The number of threads to use during read classification. By default, it is 4.
- `-rlns/--read-lines`: The number of lines in the read fasta/fastq file that are dedicated to a single read. By default, it is 2 (read name + the read itself), but this could be 4 in some cases.

### Syntax Examples

#### Example 1: MinHash-based Screen with corresponding read sets (of reads with average length 10k, error rate 1%)
To generate a MinHash-based screen for a set of genomes, and then run a corresponding set of reads (i.e. one read set for each genome) against the screen, all while using the default hash function, please run the following:

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length = 10000> -re <Expected Read Error = 0.01> -tm <Number of target matches per read>
```

#### Example 2: Minimizer-based Screen with corresponding read sets, with a specified number of reads loaded at a time, and the mmh3 hash function

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -s M -c -cs <Number of reads to be loaded at a time> -hf mmh3
```

#### Example 3: Uniform-sampling Screen of a fixed size with corresponding read sets, with a specified number of reads loaded at a time, the default hash function, with read logging active

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -f <Screen Size> -c -cs <Number of reads to be loaded at a time> -rlg -rlc <Location to save the read logs>
```

#### Example 4: MinHash-based Screen with novel/unmatched read sets, with a specified number of reads loaded at a time, the default hash function, with read logging active

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -um  -c -cs <Number of reads to be loaded at a time> -rlg -rlc <Location to save the read logs>
```

#### Example 5: Minimizer-based Screen, without any reads to be classified, using the 128-bit mmh3 hash function
```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -hf mmh3_128 -so -sl <Location to save the generated screens>
```

### Further implementation details

#### Hash Function Options

At this point, the default hash function is Java's built in `hashCode()`. For the following section, if no hash function is specified, then this will be used. However, to specify a particular function, please use the following parameters:

- `h`: Java's `hashCode()`. Default option.
- `mmh3`: Google Guava's 32-bit `MurmurHash3` implementation.
- `mmh3_128`: Google Guava's 128-bit `MurmurHash3` implementation.

The selected hash function will be used to generate the screen and to classify the reads.

I will be adding more hash functions soon!

#### Threads

The number of threads determines the number of reads processed at once. Going forward, I hope to parallelize the process of screen generation too.

#### Read Loading

As specified above, if `-c` is used, then reads will be loaded in chunks of the size specified by `-cs`. This is recommended for large read files, so that too much memory is not used.

#### Screen/Sketch Generation Only

Currently, there is only an option to save the generated screens. I am working on adding an option to load and use saved screens.

#### Read Logging

When running in matched mode (and the `-rlg` option is used), read logs will be saved to the location specified by `-rlc`. The logs will be a single line, with four numbers - the true source of the read, the predicted source, and the matches with each. The log files will be named `<ReadSet>_<ReadNumber>.log`.

When running in unmatched mode (`-um`),read logs are generated by default, and saved to the location specified by `-rlc`. The logs will be a single line, with each number corresponding the number of matches between the current read and each element of the screen. The log files will be named `<ReadSet>_<ReadNumber>.log`.

Python helper scripts for aggregating these log files and summarizing results are included, and are described below.

#### Using External JARs

As I am not using Maven, I have manually included the JARs this project will use. The foremost of these is Google's `Guava`, which gives me access to hash function implementations, optimized data structures, and other nifty features. Other include JARs are:
- `commons.cli`: This is the library used for command line argument parsing.

As I add more to this folder, I will update the README to include descriptions of each of them!

# Analyzing Results

This is still new, and I will update this over the coming weeks.

In the `analysis` folder, you can find some scripts that are useful for analyzing the results of the screening process.
- `extract_results.py`: Use `python3 extract_results.py <Path to output file>` to generate a summary of the experiment. This will give you the total number of correctly and incorrectly classified reads, as well as a few other metrics. It will also generate a histogram with the classification accuracy of reads from each organism.
- `aggregate_read_logs.py`: This is for analyzing results at the end of a matched run. Use `python3 aggregate_read_logs.py <Path to folder with read logs> <numbers of members in the community>` to get a matrix showing how many reads were correctly or incorrectly classified to each member of the community. The matrix is interpreted as follows - the value at` Matrix[x][y]` is the number of reads from organism `x` that were classified as being from organism `y`. Note that this can only be used in matched mode, if read logging was enabled.
- `aggregate_classification_logs.py`: This is for analyzing results at the end of an unmatched run. Use `python3 aggregate_classification_logs.py <Path to folder with read logs>` to generate a series of histograms showing the breakdown of how many reads were classified to each genome in the screen for each readset that was classified. The folder containing these log files will be the location stored in the `READ_LOCATION` parameter in `Settings.java`.
