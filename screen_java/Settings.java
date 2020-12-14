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

  // IO
  static String GENOME_FOLDER;
  static String READS_FOLDER;
  static String OUTPUT_FILE;

  // GENOME AND READS
  static String[] READ_SETS;
  static String[] GENOMES;

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

  // MONITOR FOR ERROR
  static boolean BAD_INPUT;

  // CHUNK BASED READ LOADING
  static boolean IN_CHUNKS;
  static int CHUNK;
  static boolean CHUNK_UPDATES;

  // Only screen generation
  static boolean SCREEN_ONLY;
  static boolean INDIVIDUAL_SCREENS;
  static String SCREEN_LOCATION;

  static void parseArgs(String[] args) throws Exception {

    // TODO: Make these parameters?

    // Experiment parameters
    K = 21;
    NUM_THREADS = 4;
    READ_LINES = 2;

    // Read Loading
    IN_CHUNKS = true;
    CHUNK = 2000;
    CHUNK_UPDATES = false;

    // Defaults
    BAD_INPUT = false;
    FIXED = false;
    FIXED_SIZE = 0;

    // Screen only option
    SCREEN_ONLY = false;
    INDIVIDUAL_SCREENS = true;

    // Read logging
    READ_LOGGING = false;

    // TODO - add a load screen from file option

    // Screen type
    FIXED = false;
    MINHASH = false;
    UNIFORM = false;
    MINIMIZER = false;

    // For debugging - TODO: maybe remove
    TRACK_MINIMIZERS = false;
    TRACK_READS = false;
    TRACK_MISCLASSIFIED = false;

    // Hash options
    // TODO - Add more, and update "getHashName" in ScreenGenerator
    List<String> hashTypes = Arrays.asList("h", "mmh3", "mmh3_128");
    // Default is built-in hashcode
    HASH_TYPE = "h";

    // Number of arguments
    int argNum = args.length;

    // Too few
    if (argNum < 6) {
      BAD_INPUT = true;
      return;
    }

    // If Screen-Only - TODO: add option for fixed size screens
    if (args[0].equals("S") || args[0].equals("s")) {

      SCREEN_ONLY = true;
      GENOME_FOLDER = args[1];
      SCREEN_LOCATION = args[2];
      READ_LENGTH = Integer.parseInt(args[3]);
      READ_ERROR = Double.parseDouble(args[4]);

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

      TARGET_MATCHES = Integer.parseInt(args[5]);

      if (argNum == 6) {
        MINHASH = false;
      } else {
        if (args[6].equals("m")) {
          MINIMIZER = true;
        } else if (args[6].equals("u")) {
          UNIFORM = true;
        }
      }

    // Default option
    } else {

      GENOME_FOLDER = args[0];
      OUTPUT_FILE = args[1];
      READS_FOLDER = args[2];
      READ_LENGTH = Integer.parseInt(args[3]);
      READ_ERROR = Double.parseDouble(args[4]);

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

      // PARSE TO FIGURE OUT WHICH SCREEN TO USE

      // Number of target Matches - initialized to zero, then read in
      TARGET_MATCHES = 0;

      // Absolute base option - MinHash Screen
      if (argNum == 6){
        TARGET_MATCHES = Integer.parseInt(args[5]);
        MINHASH = true;

      // MinHash Screen with specified hash function
      } else if ((argNum == 7) && hashTypes.contains(args[6])){
        TARGET_MATCHES = Integer.parseInt(args[5]);
        HASH_TYPE = args[6];
        MINHASH = true;

      // Remaining options - 7 parameters, 8 if specified hash function
      } else {

        // Get the hash function to be used
        if (argNum == 8 && hashTypes.contains(args[7])){
          HASH_TYPE = args[7];
        } else if (argNum == 8 && !hashTypes.contains(args[7])){
          BAD_INPUT = true;
          System.out.println("Invalid hash function specified - using default Java hash.");
        }

        // Uniformly Sampled Screen
        if (args[6].equals("u")){
          TARGET_MATCHES = Integer.parseInt(args[5]);
          UNIFORM = true;

        // Minimizer Screen, with calculated window size
        } else if (args[6].equals("m")) {
          TARGET_MATCHES = Integer.parseInt(args[5]);
          MINIMIZER = true;

        // Minimizer Screen, with provided window size
        } else if (args[5].equals("m")) {
          FIXED = true;
          // Window size
          FIXED_SIZE = Integer.parseInt(args[6]);
          MINIMIZER = true;

        // Fixed Size MinHash Screen
        } else if (args[5].equals("f")) {
          // Screen size
          FIXED = true;
          FIXED_SIZE = Integer.parseInt(args[6]);
          MINHASH = true;

        // Invalid parameters
        } else {
          BAD_INPUT = true;
          System.out.println("Invalid input parameters - please read the README!");
        }
      }
    }

  }
}
