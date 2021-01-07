import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeSet;
import java.util.*;
import java.lang.*;
import java.util.stream.Collectors;
import java.util.Collection;
import java.nio.charset.StandardCharsets;
import org.apache.commons.cli.*;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

// TODO - fix this eventually.
public class Settings {

  // Key Parameters
  static int K;
  // Number of lines each read takes
  static int READ_LINES;

  // READ PARAMS
  static int READ_LENGTH;
  static double READ_ERROR;
  static int TARGET_MATCHES;

  // CLASSIFICATION PARAMETERS
  static int THRESHOLD;
  static int NUM_THREADS;
  static int FIXED_SIZE;
  static String HASH_TYPE;

  // I/O
  static String GENOME_FOLDER;
  static String READS_FOLDER;
  static String OUTPUT_FILE;

  // GENOME AND READS
  static String[] READ_SETS;
  static String[] GENOMES;

  // KNOWN CLASSIFICATION VS NOVEL CLASSIFICATION
  static boolean MATCHED_READS_GENOMES;

  // CHOICE OF SCREEN
  static boolean FIXED;
  static boolean MINHASH;
  static boolean UNIFORM;
  static boolean MINIMIZER;

  // DEBUGGING
  static boolean TRACK_MINIMIZERS;
  static boolean TRACK_READS;
  static boolean TRACK_MISCLASSIFIED;

  // SAVE READ RESULTS
  static boolean READ_LOGGING;
  static String READ_LOCATION;

  // MONITOR FOR ERROR
  static boolean BAD_INPUT;

  // CHUNK BASED READ LOADING
  static boolean IN_CHUNKS;
  static int CHUNK;
  static boolean CHUNK_UPDATES;

  // SCREEN GENERATION WITHOUT READ SCREENING
  static boolean SCREEN_ONLY;
  static boolean INDIVIDUAL_SCREENS;
  static String SCREEN_LOCATION;

  static void parseArgs(String[] args) throws Exception {

    // Input options
    Options options = new Options();

    Option k = new Option("k", "kmer", true, "K-mer Size (default = 21)");
    k.setRequired(false);
    options.addOption(k);

    Option num_threads = new Option("nt", "num-threads", true, "Number of threads to use (default = 4)");
    num_threads.setRequired(false);
    options.addOption(num_threads);

    Option read_lines = new Option("rlns", "read-lines", true, "Number of lines per read (default = 2)");
    read_lines.setRequired(false);
    options.addOption(read_lines);

    Option in_chunks = new Option("c", "chunks", false, "Load read in chunks (default = false)");
    in_chunks.setRequired(false);
    options.addOption(in_chunks);

    Option chunk_size = new Option("cs", "chunk-size", true, "Number of reads to load at once (default = 2000)");
    chunk_size.setRequired(false);
    options.addOption(chunk_size);

    Option chunk_updates = new Option("cu", "chunk-updates", false, "Whether to print results after each chunk (default = false)");
    chunk_updates.setRequired(false);
    options.addOption(chunk_updates);

    Option screen_only = new Option("so", "screen-only", false, "Only generate screen (default = false)");
    screen_only.setRequired(false);
    options.addOption(screen_only);

    Option individual_screens = new Option("is", "individual-screens", false, "If only doing screen generation save each genome separately (default = true)");
    individual_screens.setRequired(false);
    options.addOption(individual_screens);

    Option screen_location = new Option("sl", "screen-location", true, "Where to save screen (default = ./screens)");
    screen_location.setRequired(false);
    options.addOption(screen_location);

    Option read_logging = new Option("rlg", "read-logging", false, "Save read logs (default = true)");
    read_logging.setRequired(false);
    options.addOption(read_logging);

    Option read_location = new Option("rlc", "read-location", true, "Directory for read logs (default = ./logs)");
    read_location.setRequired(false);
    options.addOption(read_location);

    // Main parameters
    Option genome_folder = new Option("g", "genome", true, "Directory containing genomes");
    genome_folder.setRequired(true);
    options.addOption(genome_folder);

    Option read_folder = new Option("r", "reads", true, "Directory containing reads");
    read_folder.setRequired(true);
    options.addOption(read_folder);

    Option read_length = new Option("rl", "read-length", true, "Read Lengths");
    read_length.setRequired(false);
    options.addOption(read_length);

    Option read_error = new Option("re", "read-error", true, "Read Error");
    read_error.setRequired(false);
    options.addOption(read_error);

    Option target_matches = new Option("tm", "target-matches", true, "Read Error");
    target_matches.setRequired(false);
    options.addOption(target_matches);

    Option screen_type = new Option("s", "screen-type", true, "Which screen-generation approach to use (default = MinHash, options = minhash, [u]niform, [m]inimizer)");
    screen_type.setRequired(false);
    options.addOption(screen_type);

    Option output_file = new Option("o", "output", true, "Where to save experiment output (default = ./)");
    output_file.setRequired(false);
    options.addOption(output_file);

    Option hash_type = new Option("h", "hash", true, "Which hash function to use (default = Java's built in Hashcode, options = hashcode, [mmh3], [mmh3-128])");
    hash_type.setRequired(false);
    options.addOption(hash_type);

    Option fixed = new Option("f", "fixed", true, "If the screen size should be fixed, then specify what size it the screens should be");
    fixed.setRequired(false);
    options.addOption(fixed);

    Option unmatched = new Option("um", "unmatched", false, "Genomes and reads do not match (default = genomes and reads match)");
    unmatched.setRequired(false);
    options.addOption(unmatched);

    // Parse

    CommandLineParser parser = new DefaultParser();
    HelpFormatter formatter = new HelpFormatter();
    CommandLine cmd;

    try {
        cmd = parser.parse(options, args);
    } catch (ParseException e) {
        System.out.println(e.getMessage());
        formatter.printHelp("Sketching", options);
        System.exit(1);
        return;
    }

    // Tracking parameters
    BAD_INPUT = false;

    // Experiment-wide parameters
    if (cmd.hasOption("k")) {
      K = Integer.parseInt(cmd.getOptionValue("k"));
    } else {
      K = 21;
    }

    if (cmd.hasOption("nt")) {
      NUM_THREADS = Integer.parseInt(cmd.getOptionValue("nt"));
    } else {
      NUM_THREADS = 4;
    }

    if (cmd.hasOption("rlns")) {
      READ_LINES = Integer.parseInt(cmd.getOptionValue("rlns"));
    } else {
      READ_LINES = 2;
    }

    // Chunk loading options
    if (cmd.hasOption("c")) {
      IN_CHUNKS = true;
    } else {
      IN_CHUNKS = false;
    }

    if (cmd.hasOption("cs")) {
      CHUNK = Integer.parseInt(cmd.getOptionValue("cs"));
    } else {
      CHUNK = 2000;
    }

    if (cmd.hasOption("cu")) {
      CHUNK_UPDATES = true;
    } else {
      CHUNK_UPDATES = false;
    }

    // Screen-generation only options
    if (cmd.hasOption("so")) {
      SCREEN_ONLY = true;
    } else {
      SCREEN_ONLY = false;
    }

    if (cmd.hasOption("is")) {
      INDIVIDUAL_SCREENS = true;
    } else {
      INDIVIDUAL_SCREENS = false;
    }

    if (cmd.hasOption("sl")) {
      SCREEN_LOCATION = cmd.getOptionValue("sl");
    } else {
      SCREEN_LOCATION = "./screens";
    }

    // Read Logging Options
    if (cmd.hasOption("rlg")) {
      READ_LOGGING = true;
    } else {
      READ_LOGGING = false;
    }

    if (cmd.hasOption("rlc")) {
      READ_LOCATION = cmd.getOptionValue("rlc");
    } else {
      READ_LOCATION = "./logs";
    }

    // KEY PARAMETERS

    // Folder containing genomes
    GENOME_FOLDER = cmd.getOptionValue("g");

    // Folder containing reads
    READS_FOLDER = cmd.getOptionValue("r");

    // Output File
    OUTPUT_FILE = cmd.getOptionValue("o");

    // Screen type
    FIXED = false;
    MINHASH = false;
    UNIFORM = false;
    MINIMIZER = false;

    String st = cmd.getOptionValue("s");
    if (st == null){
      MINHASH = true;
    } else if (st.equals("m") || st.equals("M")) {
      MINIMIZER = true;
    } else if (st.equals("u") || st.equals("U")) {
      UNIFORM = true;
    } else {
      System.out.println("Invalid Screen Generation approach specified - using MinHash for screen generation");
      MINHASH = true;
    }

    if (cmd.hasOption("f")) {
      FIXED = true;
      FIXED_SIZE = Integer.parseInt(cmd.getOptionValue("f"));
    } else {
      FIXED = false;
      FIXED_SIZE = 0;
    }

    // Read attributes
    if (cmd.hasOption("rl")) {
      READ_LENGTH = Integer.parseInt(cmd.getOptionValue("rl"));
    } else {
      if (FIXED) {
        READ_LENGTH = 0;
      } else {
        System.out.println("Invalid Read Length (-rl/--read-length) specified.");
        BAD_INPUT = true;
        READ_LENGTH = 10000;
      }
    }

    if (cmd.hasOption("re")) {
      READ_ERROR = Double.parseDouble(cmd.getOptionValue("re"));
    } else {
      if (FIXED) {
        READ_ERROR = 0;
      } else {
        System.out.println("Invalid Read Error (-re/--read-error) specified.");
        BAD_INPUT = true;
        READ_ERROR = 0.01;
      }
    }

    // Read attributes
    if (cmd.hasOption("tm")) {
      TARGET_MATCHES = Integer.parseInt(cmd.getOptionValue("tm"));
    } else {
      if (FIXED) {
        TARGET_MATCHES = 0;
      } else {
        System.out.println("Invalid TargetMatches (-tm/--target-matches) specified.");
        BAD_INPUT = true;
        TARGET_MATCHES = 30;
      }
    }

    // Hash Function being used
    List<String> hashTypes = Arrays.asList("h", "mmh3", "mmh3_128");
    if (cmd.hasOption("h")) {
      String selected_hash = cmd.getOptionValue("h");
      if (hashTypes.contains(selected_hash)){
        HASH_TYPE = selected_hash;
      } else {
        System.out.println("Invalid hash function specified - using default Java hash.");
        HASH_TYPE = "h";
      }
    } else {
      HASH_TYPE = "h";
    }

    if (cmd.hasOption("um")) {
      MATCHED_READS_GENOMES = false;
    } else {
      MATCHED_READS_GENOMES = true;
    }

    // Debugging Parameters
    TRACK_MINIMIZERS = false;
    TRACK_READS = false;
    TRACK_MISCLASSIFIED = false;

    // Check if anything is broken by now
    if (BAD_INPUT) {
      return;
    }

    // Screen Generation only - TODO: Add option for fixed size screens
    if (SCREEN_ONLY) {
      // Get a sorted list list of all genomes in the input folder
      List<String> genomes_list = new ArrayList<String>();
      File dir = new File(GENOME_FOLDER);
      for (File file : dir.listFiles()) {
        if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
          genomes_list.add(file.getName());
        }
      }
      Collections.sort(genomes_list);
      GENOMES = new String[genomes_list.size()];
      GENOMES = genomes_list.toArray(GENOMES);

    // Default option - screen generation and read screening
    } else {

      // Get a sorted list list of all genomes in the input folder
      List<String> genomes_list = new ArrayList<String>();
      File dir = new File(GENOME_FOLDER);
      for (File file : dir.listFiles()) {
        if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
          genomes_list.add(file.getName());
        }
      }
      Collections.sort(genomes_list);
      GENOMES = new String[genomes_list.size()];
      GENOMES = genomes_list.toArray(GENOMES);

      // Get a sorted list of all readsets in the input folder
      List<String> reads_list = new ArrayList<String>();
      File readDir = new File(READS_FOLDER);
      for (File readFile : readDir.listFiles()) {
        if (readFile.getName().endsWith((".fasta")) || readFile.getName().endsWith((".fna")) || readFile.getName().endsWith((".fa"))) {
          reads_list.add(readFile.getName());
        }
      }
      Collections.sort(reads_list);
      READ_SETS = new String[reads_list.size()];
      READ_SETS = reads_list.toArray(READ_SETS);
    }

  }
}
