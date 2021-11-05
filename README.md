# Overview
All the code used for "Analyzing sketching and sampling approaches for fast and accurate long read classification".

This README contains details about the code and how to use it, information about the auxiliary scripts used to analyze the results, and then details about how to generate simulated reads and access the data we used. This repository

I will continue to add to it and update the README to reflect any changes.

# Sections
- [Pre-prints, Posters and Presentations](#pre-prints--posters-and-presentations)
- [Implementation Details](#implementation-details)
    + [Brief overview of the approach](#brief-overview-of-the-approach)
    + [Running the code and explanation of parameters](#running-the-code-and-explanation-of-parameters)
    + [Syntax Examples](#syntax-examples)
    + [Further implementation details](#further-implementation-details)
      - [Hash Function Options](#hash-function-options)
      - [Threads](#threads)
      - [Specifying input files](#specifying-input-files)
      - [Read Loading](#read-loading)
      - [Screen/Sketch Generation Only](#screen-sketch-generation-only)
      - [Read Logging](#read-logging)
        * [Matched Mode](#matched-mode)
        * [Unmatched Mode](#unmatched-mode)
      - [Using External JARs](#using-external-jars)
- [Analyzing Results](#analyzing-results)
- [Generating Reads](#generating-reads)
    + [Read Length Options (Not recommended)](#read-length-options--not-recommended-)
- [Data availability](#data-availability)

# Pre-prints, Posters and Presentations

- The pre-print is out! You can find here, on bioRxiv: https://www.biorxiv.org/content/10.1101/2021.11.04.467374v1

- [The poster of this work at Genome Informatics 2021](https://drive.google.com/file/d/1CUPuw4Wv_yJiZ4jZkDbVVVF8dx3MMoC6/view?usp=sharing). Link to the lightning talk for this will be posted soon!

- [The poster of this work at RECOMB 2021](https://drive.google.com/file/d/1rNIQLJqiZzNdNJfhrTTExDxzFXM7-LLD/view?usp=sharing) (where it won the award for "Best Poster"!), and [my poster presentation](https://drive.google.com/file/d/1bjxqqSY75X8wEqR6_EUYMB_PbVK78IAL/view?usp=sharing).
- [The lightning talk for this poster from RECOMB 2021](https://youtu.be/MMX2kbkRmYI).

- [The poster of this work at Biological Data Science 2020](https://drive.google.com/file/d/1gFN_F26f3UZwSZawbVyWXuP3DnGB9_sH/view?usp=sharing).

# Implementation Details
The tool is built in java, and provides a way to generate a "screen" (sketched representation) of an input set of genomes, and either save the screen or classify input reads against this screen. For read classification, these may be reads drawn from the same set of genomes in the screen (either simulated using the approach above, or from another source), or a totally different set of reads.

There are three distinct approaches for screen-generation: a MinHash-based approach (with a few variants), a Minimizer-based approach, and a Uniform sampling approach. For each of these approaches, the size of the screen can either be calculated based on user-specified input parameters (the expected read length, the expected read error, the sizes of the genomes, and the number of shared samples each read should aim to have with its correct genome), or can be specified ahead of time by the user.

To see details of how to do any of the experiments outlined above, or how to adjust certain options/settings, please read the later sections of the README.

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
- `-g/--genome <Directory>`: The directory containing the genomes that are to be sketched into a screen. Only not necessary if we are using pre-generated screens (`-ls/--load-screens`). In this situation, the pre-generated screens stored at `-sl` will be used. You can also pass in a file containing a list of files (with absolute paths) instead of specifying a folder - for more information on this, see [this section](#specifying-input-files).
- `-r/--reads <Directory>`: The directory containing the reads that are to be classified against the screen. Only not necessary if we are in "screen-only" mode (`-so/--screen-only`). You can also pass in a file containing a list of files (with absolute paths) instead of specifying a folder - for more information on this, see [this section](#specifying-input-files).
- `-rl/--read-length <Integer>`: The expected read length of reads that will be classified against this screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-re/--read-error <Double>`: The expected error rate of reads that will be classified against this screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-tm/--target-matches <Integer>`: The target number of matches between a read and its correct source in the screen. Not necessary if a fixed size screen (`-f`) is being generated.
- `-s/--screen-type`: The method of screen-generation to be used for this screen. Choose between MinHash-based (default), [m]inimizer based, [u]niform sampling or [e]xhaustive.
- `-hf/--hash <Hash Type>`: The hash function to use throughout the process. Choose between Java's built in hashcode (default), [mmh3] (MurmurHash3) or [mmh3_128] (the 128-bit variant of MurmurHash3).
- `-th/--threshold <Integer>`: This is the number of matches a read must have with a source genome to be classified. By default, it is 5 - however, this should be adjusted based on the level of confidence we want in classification calls. Reads below this threshold will be flagged as insufficient. This value should always be below the number of target matches, and usually just a fraction of the target matches (so that reads are classified as often as possible).

Matched Reads/Genomes vs Classification without a ground truth:
- `-um/--unmatched`: Use this parameter if the reads and the genomes do not correspond - by default, the code will assume that there is a matching read set for each genome (generated using the simulator above), and will proceed with classification assuming there is some ground truth. If this option is selected, detailed classification logs will be saved for each read, and no accuracy metrics will be printed. NOTE: Over time, this will become the default, and a parameter will be needed to indicate the genomes and reads are matched. In this mode, read logging is enabled by default.

Additional Options for MinHash:
NOTE: These sections are still under development.
- `-wmh/--weighted-minhash`: Use this flag to enable weighted minhash. By default, it is false.
- `-wmhw/--weighted-minhash`: If `wmh` is used, this can be used to specify the multiplier given to unique k-mers (i.e. k-mers that only occur in one genome in the screen). By default it is one (as all k-mer weights are the total number of genomes - the number of genomes they occur in), but with this option unique k-mers can be weighted higher, helping to break ties.
- `-omh/--order-minhash`: Use this flag to enable order minhash. By default, it is false.
- `-omhl/--order-minhash-len`: If `-omh` is used, then this will specify the number of k-mers whose order should be retained relative to each other. By default, this is 3.

Fixed Size Screens:
- `-f/--fixed <Integer>`: Use if the screen size should not be calculated, but instead the specified size/window size should be used for all screens. For uniform and MinHash screens, this will be the sketch size. For Minimizer-based screens, this will be the window size. This is not available for exhaustive screens, which just use all hashes in a given genome.

Read Loading:
- `-c/--chunks`: Use if a specified number of reads should loaded at a time (i.e. in "chunks"), instead of one file at a time.
- `-cs/--chunk-size <Integer>`: If `-c` is used, then `-cs` will specify the "chunk" size.
- `-cu/--chunk-updates`: If `-c` is used, then enabling `-cu` will print an update after each chunk is processed.

Read Logging:
- `-rlg/--read-logging`: Use if read-logging should be enabled. This is enabled by default in unmatched mode (`-um`).
- `rlc/--read-location <Directory>`: The directory that read logs should be saved to. They will be named with the naming convention `<Readset Number>_<Read_Number>.log` - for example, `3_23.log` is the log for the 24th read from the 3rd read set.

Generating Screens without doing classification:
- `-so/--screen-only`: Use if you only want to generate the screen, but not classify any reads.
- `-sl/--screen-location <Directory>`: Used for both screen-generation and loading pre-generated screens. Specifies the location where the generated screens will be saved. Screens will be saved in `.bin` files with names matching the input genome files. If weighted minhash or order minhash are being used, then an additional `.bin` file will be saved with the weights/order, and will be loaded when used. Finally, a `params.txt` file with the details of the experiment will be generated - you will need to use the information in this file when loading the saved screen.

Load Pre-Generated Screen:
- `-ls/--load-screen`: Use this if you want to use a pre-generated screen, and supply the input reads to classify against it. Please remember to set matched/unmatched (`-um`) depending on the type of classification you want to do, and please specify the error rate, read length and target matches used when generating this screen, as well as the weighted or order minhash parameters if either was used. The details for the screens that are loaded can be found in the same directory in the `params.txt` file, so refer to that to get the correct values.
- `-sl/--screen-location <Directory>`: Used for both screen-generation and loading pre-generated screens. Specifies the location where the pre-generated screens will be loaded from. Pre-generated screens will be loaded from `.bin` files with names matching the input genome files.

Save Screen while doing regular classification:
- `-ss/--save-screen`: Use this to save screens after generation, but before the classification stage. Use the `-sl` parameter to specify the location it is saved to.

Experiment Parameters
- `-k/--kmer <Integer>`: The k-mer size to use. By default, it is 21.
- `-nt/--num-threads <Integer>`: The number of threads to use during read classification. By default, it is 4.
- `-rlns/--read-lines`: The number of lines in the read fasta/fastq file that are dedicated to a single read. By default, it is 2 (read name + the read itself), but this could be 4 in some cases.

Clustering-based Approach:
NOTE: This section is still under development.
- `-ct/--cluster`: Use this flag to enable the clustering approach, where the generated screens will be clustered and arranged into a tree. By default, this is false, and regular classification will occur.
- `-css/--cluster-sketch-size`: Used to specify the size of the sketches used to cluster the genomes. By default is set to 100, but can be set to any value as long as clustering is enabled.
- `dt/--downsample-type`: Used to specify the downsampling format we want to use in cluster generation. By default, sketches are `[n]`ot downsampled, but other options include downsampling by a `[c]`onstant factor or scaled for `[h]`eight.
- `df/--downsample-factor`: If a downsampling format is specified, you can specify the factor used when downsampling. By default this is one, but any integer may be provided.

Read Filtering:
- `-fr/--filter-reads`: Use this flag to filter out low quality reads. For simulated reads, this will remove reads that have a large number of N's. In the future, this will filter out low quality reads.

### Syntax Examples

#### Example 1: MinHash-based Screen with corresponding read sets (of reads with average length 10k, error rate 1%)
To generate a MinHash-based screen for a set of genomes, and then run a corresponding set of reads (i.e. one read set for each genome) against the screen, all while using the default hash function, please run the following:

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length = 10000> -re <Expected Read Error = 0.01> -tm <Number of target matches per read>
```

#### Example 2: Minimizer-based Screen with corresponding read sets, with a specified number of reads loaded at a time, and the mmh3 hash function

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -s m -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -s M -c -cs <Number of reads to be loaded at a time> -hf mmh3
```

#### Example 3: Uniform-sampling Screen of a fixed size with corresponding read sets, with a specified number of reads loaded at a time, the default hash function, with read logging active

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -s u -f <Screen Size> -c -cs <Number of reads to be loaded at a time> -rlg -rlc <Location to save the read logs>
```

#### Example 4: MinHash-based Screen with novel/unmatched read sets, with a specified number of reads loaded at a time, the default hash function, using 4 threads

```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -um  -c -cs <Number of reads to be loaded at a time> -rlc <Location to save the read logs> -nt 4
```

#### Example 5: Minimizer-based Screen, without any reads to be classified, using the 128-bit mmh3 hash function
```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -s m -rl <Expected Read Length> -re <Expected Read Error> -tm <Number of target matches per read> -hf mmh3_128 -so -sl <Location to save the generated screens>
```

#### Example 6: Use a clustered uniform-sampling approach with corresponding read sets, with clustering done using a sketch size of 250, with a specified number of reads  loaded at a time, the default hash function, and read logging active
```
java -cp screen_java:jars/\* Main -g <Genomes Directory> -r <Reads Directory> -o <Log File location> -s u -f <Screen Size> -c -cs <Number of reads to be loaded at a time> -rlg -rlc <Location to save the read logs> -ct -css 250
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

The number of threads parameter determines the number of reads processed at once. Going forward, I hope to parallelize the process of screen generation too.

#### Specifying Input Files

Instead of passing in a folder containing the genomes or reads, you may instead pass in a file containing a list of files to be used. To generate such a file you can navigate to the folder you want, and use the command `ls -d "$PWD"/* > genomes.txt` to generate a list of files in the folder, listed with absolute path. This option allows you to pick and choose which files in a folder you want to include, instead of just using the entire folder.

#### Read Loading

As specified above, if `-c` is used, then reads will be loaded in chunks of the size specified by `-cs`. This is recommended for large read files, so that too much memory is not used.

#### Screen/Sketch Generation Only

You can generate and save a screen for a later experiment, or load an existing screen. To do these, please use the `-so` or `-ls` parameters respectively, with `-sl` specifying which location the screens are to be saved to/loaded from.

#### Read Logging

##### Matched Mode
When running in matched mode (and the `-rlg` option is used), read logs will be saved to the location specified by `-rlc`. The format of the log files varies depending on whether or not clustering is used. If clustering is not used, the logs will be a two lines:
- An array of integers, with each number corresponding the number of matches between the current read and each element of the screen,
-  The true source of the read, noted as an integer value (e.g. if it is `5`, then the 6th input genome is the source).

If clustering is used, then the log file will contain five lines:

- The number of genomes in the screen,
- The path down the tree taken to classify the read,
- The predicted source genome (the index of it),
- The number of matches with that source genome,
- The true source of the read

The log files will be named `<ReadSet>_<ReadNumber>.log`.

##### Unmatched Mode
When running in unmatched mode (`-um`),read logs are generated by default, and saved to the location specified by `-rlc`. As in the matched situation, the format of the log files will vary depending on whether or not clustering was used. Without clustering, the logs will be a single line, with each number corresponding the number of matches between the current read and each element of the screen. With clustering, the format is almost identical to that specified in the matched section above, with the only change being that there is no fifth line (as the true source is not known).

The log files will be named `<ReadSet>_<ReadNumber>.log`.

Python helper scripts for aggregating these log files and summarizing results are included, and are described below.

#### Using External JARs

As I am not using Maven, I have manually included the JARs this project will use. The foremost of these is Google's `Guava`, which gives me access to hash function implementations, optimized data structures, and other nifty features. Other include JARs are:
- `commons.cli`: This is the library used for command line argument parsing.
- `hierarchical-clustering`: This is the package used for clustering genome sketches to construct the sketch tree.
- `clust4j`: An alternate option for clustering - currently not in use.

As I add more to this folder, I will update the README to include descriptions of each of them!

# Analyzing Results

This is still new, and I will update this as things get added!

In the `analysis` folder, you can find some scripts that are useful for analyzing the results of the screening process.
- `extract_results.py`: Use `python3 extract_results.py <Path to output file>` to generate a summary of the experiment. This will give you the total number of correctly and incorrectly classified reads, as well as a few other metrics. It will also generate a histogram with the classification accuracy of reads from each organism.
- `aggregate_read_logs.py`: This is for analyzing results at the end of a matched run. Use `python3 aggregate_read_logs.py <Path to folder with read logs>` to get read set-by-read set breakdown of the classification accuracy, as well as a matrix showing how many reads were correctly or incorrectly classified to each member of the community. The matrix is interpreted as follows - the value at` Matrix[x][y]` is the number of reads from organism `x` that were classified as being from organism `y`. Note that this can only be used in matched mode, if read logging was enabled. This file will also generate a number of plots and visualizations, such as the breakdown of how each readset was classified, a visualization of the classification matrix, and a distribution of how many matches reads had with their source genome.
- `aggregate_classification_logs.py`: This is for analyzing results at the end of an unmatched run. Use `python3 aggregate_classification_logs.py <Path to folder with read logs>` to generate a series of histograms showing the breakdown of how many reads were classified to each genome in the screen for each readset that was classified. The folder containing these log files will be the location stored in the `READ_LOCATION` parameter in `Settings.java`. As with the previous file, there will also be a number of plots and visualizations generated.

# Generating Reads
The read simulator is a modified version of [Melanie Kirsche's read simulator](https://github.com/schatzlab/centroTools/tree/master/java). To generate reads for a set of genomes, run `generate_reads.sh` using the following syntax:
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

By default, the normal distribution will be used. However, to use any of the other four options, simply add the appropriate string (`XL/E/EM/EL`) as an additional parameter when generating reads - for example, running the command above with `XL` added (`./generate_reads.sh readsim MBARC_ZYMO reads_MBARC_ZYMO 0.0034 0.0033 0.0033 10000 10`) will now generate reads of exactly length 10KB.

# Data availability

ZYMO: https://www.zymoresearch.com/collections/zymobiomics-microbial-community-standards

MBARC-26: https://www.nature.com/articles/sdata201681

CGR: https://www.nature.com/articles/s41587-018-0008-8

PacBio HiFi Microbiome reads: https://www.pacb.com/blog/data-release-human-microbiome-samples-demonstrate-advances-in-hifi-enabled-metagenomic-sequencing/
