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

public class ReadScreener {
  // Key Variables
  int readLen;
  double readErr;
  int targetMatches;
  String genomeFolder;
  String readFolder;
  String[] genomeNames;
  String[] readSets;
  int numGenomes;
  ArrayList<HashSet<Integer>> sketch_hash;
  int window;

  public static void main(String[] args) throws Exception
  {
    // ----- PARAMETERS -----
    int k = 21;

    // Debugging parameters - off by default
    boolean track_minimizers = false;
    boolean track_reads = false;

    // Number of arguments
    int argNum = args.length;

    // Genome Folder
    String gf = args[0];

    // Output file
    String ofn = args[1];

    // Read Folder
    String rf = args[2];

    // Read Length
    int rl = Integer.parseInt(args[3]);

    // Read error
    double re = Double.parseDouble(args[4]);

    // Get a sorted list list of all genomes
    List<String> genomeNamesList = new ArrayList<String>();
    File dir = new File(gf);
    for (File file : dir.listFiles()) {
      if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
        genomeNamesList.add(file.getName());
      }
    }
    Collections.sort(genomeNamesList);
    String[] g = new String[genomeNamesList.size()];
    g = genomeNamesList.toArray(g);

    // Get a sorted list of all readsets
    List<String> readSetList = new ArrayList<String>();
    File readDir = new File(rf);
    for (File readFile : readDir.listFiles()) {
      if (readFile.getName().endsWith((".fasta")) || readFile.getName().endsWith((".fna")) || readFile.getName().endsWith((".fa"))) {
        readSetList.add(readFile.getName());
      }
    }
    Collections.sort(readSetList);
    String[] r = new String[readSetList.size()];
    r = readSetList.toArray(r);

    // ----- SCREENING AND CLASSIFICATION -----
    // Number of target Matches - initialized to zero, then read in
    int tm = 0;

    // Hash options
    // TODO - Add more, and update "getHashName" in ScreenGenerator
    List<String> hashTypes = Arrays.asList("h", "mmh3");
    // Default is built-in hashcode
    String hashType = "h";

    // Absolute base option - MinHash Screen
    if (argNum == 6){
      tm = Integer.parseInt(args[5]);
      MinHashScreen ms = new MinHashScreen(gf, g, rl, re, tm, k, hashType);
      ReadScreener rs = new ReadScreener(ofn, ms, gf, g, r, rf, rl, re, k, track_reads, hashType);

    // MinHash Screen with specified hash function
    } else if ((argNum == 7) && hashTypes.contains(args[6])){
      tm = Integer.parseInt(args[5]);
      hashType = args[6];
      MinHashScreen ms = new MinHashScreen(gf, g, rl, re, tm, k, hashType);
      ReadScreener rs = new ReadScreener(ofn, ms, gf, g, r, rf, rl, re, k, track_reads, hashType);


    // Remaining options - 7 parameters, 8 if specified hash function
    } else {

      // Get the hash function to be used
      if (argNum == 8 && hashTypes.contains(args[7])){
        hashType = args[7];
      } else if (argNum == 8 && !hashTypes.contains(args[7])){
        System.out.println("Invalid hash function specified - using default Java hash.");
      }

      // Uniformly Sampled Screen
      if (args[6].equals("u")){
        tm = Integer.parseInt(args[5]);
        UniformScreen us = new UniformScreen(gf, g, rl, re, tm, k, args[6], hashType);
        ReadScreener rs = new ReadScreener(ofn, us, gf, g, r, rf, rl, re, k, track_reads, hashType);

      // Minimizer Screen, with calculated window size
      } else if (args[6].equals("m")) {
        tm = Integer.parseInt(args[5]);
        MinimizerScreen mzs = new MinimizerScreen(gf, g, rl, re, tm, k, args[6], hashType);
        ReadScreener rs = new ReadScreener(ofn, mzs, gf, g, r, rf, rl, re, k, track_reads, hashType);

      // Minimizer Screen, with provided window size
      } else if (args[5].equals("m")) {
        // Window size
        int ws = Integer.parseInt(args[6]);
        MinimizerScreen mzs = new MinimizerScreen(gf, g, rl, re, k, args[5], ws, hashType);
        ReadScreener rs = new ReadScreener(ofn, mzs, gf, g, r, rf, rl, re, k, track_reads, hashType);

      // Fixed Size MinHash Screen
      } else if (args[5].equals("f")) {
        // Screen size
        int ss = Integer.parseInt(args[6]);
        MinHashScreen ms = new MinHashScreen(gf, g, rl, re, k, args[5], ss, hashType);
        ReadScreener rs = new ReadScreener(ofn, ms, gf, g, r, rf, rl, re, k, track_reads, hashType);
      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
    System.out.println("Done!");
  }

  // Main function for screening reads
  ReadScreener(String filename, ScreenGenerator sg, String gf, String[] g, String[] r, String rf, int readLength, double readError, int k, boolean track_reads, String hashType) throws Exception
  {
    int threshold = 2;

    // Store all values
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.readFolder = rf;
    this.genomeNames = g;
    this.readSets = r;
    this.sketch_hash = sg.sketch_hash;
    this.numGenomes = g.length;
    this.window = sg.window;

    // Print Sketch Sizes
    for (int a = 0; a < numGenomes; a++)
    {
      System.out.println(genomeNames[a] + " " + sketch_hash.get(a).size());
    }

    System.out.println("Screening Reads...");

    // Initialize count arrays to track accuracy
    int[] totalReads = new int[numGenomes];
    int[] correctCounts = new int[numGenomes];
    int[] misCounts = new int[numGenomes];
    int[] insufCounts = new int[numGenomes];
    int[] tieCounts = new int[numGenomes];

    // Matrix to track misclassifications
    // M[a][b] will track the number of times a is misclassified as b
    int[][] miscount_matrix = new int[numGenomes][numGenomes];

     // For each organism/element of the screen
    for (int rs = 0; rs < numGenomes; rs++)
    {
      // Initialize counts for this readset
      correctCounts[rs] = 0;
      misCounts[rs] = 0;
      insufCounts[rs] = 0;
      totalReads[rs] = 0;
      tieCounts[rs] = 0;

      // Version for automatic loading
      ArrayList<String> reads = getReads(readFolder + readSets[rs]);

      // For each read in the readset
      int numReads = reads.size();

      for (int rd = 0; rd < numReads; rd++)
      {
        // Get current read
        String currRead = reads.get(rd);
        totalReads[rs]++;

        // Get read k-mers
        ArrayList<Integer> readHashes = new ArrayList<Integer>();
        if (this.window > 0) {
          readHashes = getAllMinimizers(currRead, this.window, k, hashType);
        } else {
          readHashes = getReadKmersHash(currRead, hashType);
        }

        // Get the number of matches for this read with each genome
        int[] readScores;
        readScores = screenReadHash(sketch_hash, readHashes);

        // Get predicted source
        int pred = getMaxIndex(readScores);

        // Update counts
        if (pred == rs && readScores[pred] > threshold) {
          // Correctly classified
          correctCounts[rs]++;
          if (sameCounts(readScores, pred) > 0) {
            // Tie was broken correctly
            tieCounts[rs]++;
          }
        } else {
          if (readScores[pred] == 0){
            // Not enough to classify - should not happen
            insufCounts[rs]++;
          } else{
            // Misclassified - update counts
            misCounts[rs]++;
            miscount_matrix[rs][pred]++;
            if (readScores[pred] == readScores[rs]){
              // Tie, but this time broken incorrectly
              tieCounts[rs]++;
            }
          }
        }
      }
      // Print summary
      System.out.println(readSets[rs]);
      System.out.println(totalReads[rs] + " " + correctCounts[rs] + " " + misCounts[rs] + " " +  insufCounts[rs] + " " + tieCounts[rs]);
    }

    // Whole experiment summary
    System.out.println("Printing Results...");
    System.out.println("Number of Reads:");
    System.out.println(Arrays.toString(totalReads));
    System.out.println("Number Correctly Classified:");
    System.out.println(Arrays.toString(correctCounts));
    System.out.println("Number Incorrectly Classified:");
    System.out.println(Arrays.toString(misCounts));
    System.out.println("Number Not Classified:");
    System.out.println(Arrays.toString(insufCounts));

    // Print misclassification matrix
    printMatrix(miscount_matrix);

    //Save results
    saveResults(filename, readSets, genomeNames, sketch_hash, totalReads, correctCounts, misCounts, insufCounts, tieCounts);
  }

  // ----- HELPER FUNCTIONS FOR SCREENING ------
  // Break Read into Kmers - HASH VERSION
  ArrayList<Integer> getReadKmersHash(String read, String hashType)
  {
    int readLen = read.length();
    ArrayList<Integer> readMers = new ArrayList<Integer>();

    for (int i = 0; i < readLen-21; i++)
    {
      // Slice the kmer, and chose the canonical ones
      String readMer = read.substring(i, i+21);
      String selectedMer = getCanonical(readMer);
      readMers.add(getHash(selectedMer, hashType));
    }
    return readMers;
  }

  // Screen read kmers against screen - HASH VERSION
  int[] screenReadHash(ArrayList<HashSet<Integer>> sketch, ArrayList<Integer> readMers){

    // Number of genomes in screen
    int num = sketch.size();
    // Track number of matches
    int[] scores = new int[num];
    // For each sketch in the screen
    for (int i = 0; i < num; i++)
    {
      // Get and store overlap
      int score = getOverlapHash(sketch.get(i), readMers);
      scores[i] = score;
    }
    return scores;
  }

  // Finds number of shared elements between a sketch and a list of query kmers - HASH VERSION
  int getOverlapHash(HashSet<Integer> sketch_set, ArrayList<Integer> readMers_list)
  {
    int overlap = 0;

    // DEBUGGING - tracking number of duplicates
    int dupe = 0;
    // DEBUGGING - track the number of duplicates k-mers
    Set<Integer> already_counted  = new HashSet<Integer>();

    // Iterate through all the k-mers in the read, compare them to the sketch set
    for (int i = 0; i < readMers_list.size(); i++) {
      if (sketch_set.contains(readMers_list.get(i))) {
        if (already_counted.contains(readMers_list.get(i))) {
          // DEBUGGING - keeping track of duplicate k-mers
          dupe++;
        } else {
          already_counted.add(readMers_list.get(i));
          overlap++;
        }
      }
    }
    return overlap;
  }

  // Check if there are any ties for top value in the list
  int sameCounts(int[] matches, int topMatchIndex){
    int m = 0;
    int topMatchVal = matches[topMatchIndex];
    for (int r = 0; r < matches.length; r++){
      if (matches[r] == topMatchVal && r != topMatchIndex){
        m++;
      }
    }
    return m;
  }

  // ----- UTILITY FUNCTIONS ------
  // Essentially argmax, but with default value 0
  // Not using argmax, so that I can add a tie-breaking mechanism in there
  int getMaxIndex(int[] numArray){
    int max = 0;
    int max_pos = 0;
    for (int i = 0; i < numArray.length; i++)
    {
      if (numArray[i] > max) {
        max = numArray[i];
        max_pos = i;
      }
    }
    return max_pos;
  }

  // Get the reverse complement of a sequence
  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  // Get the lexicographically smaller of a k-mer and its reverse complement
  String getCanonical(String seq){
    String reversed_mer = reverseComplement(seq);
    String selected_mer = "";
    if (seq.compareTo(reversed_mer) > 0){
      selected_mer = reversed_mer;
    } else {
      selected_mer = seq;
    }
    return selected_mer;
  }

  // Given a string and the hash function to be used, returns the hashed sequence
  int getHash(String seq, String hashType){
    if (hashType.equals("h")){
      return seq.hashCode();
    } else {
      // TODO - fix this
      int hashVal = Hashing.murmur3_32().hashString(seq, StandardCharsets.UTF_8).asInt();
      return hashVal;
    }
  }

  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Gets all minimizers from a given string with specified window size and k
  ArrayList<Integer> getAllMinimizers(String g, int window_size, int k, String hashType){
    // Get number of kmers and windows in this string
    int num_mers = g.length() - k + 1;
    int num_windows = g.length() - window_size + 1;

    // Initialize data structures to store kmers and minimizers
    int[] string_hashes = new int[num_mers];
    HashSet<Integer> minimizers = new HashSet<Integer>();

    // Get all the k-mer hashes from the input string
    for (int i = 0; i < num_mers; i++){
      string_hashes[i] = getHash(getCanonical(g.substring(i, i+k)), hashType);
    }

    //iterate through the windows
    for (int j = 0; j < num_windows; j++){
      // Find the minimal kmer hash for the window that starts at this point
      int min_val = string_hashes[j];
      // Look through all kmers in this window
      for (int x = j + 1; x < (j + window_size - k); x++){
        int curr = string_hashes[x];
        if (curr < min_val) {
          min_val = curr;
        }
      }
      // Add this window's minimizer
      minimizers.add(min_val);
    }
    ArrayList<Integer> minimizer_list = new ArrayList<Integer>(minimizers);
    return minimizer_list;
  }

  // ----- I/O HELPER FUNCTIONS ------
  // Load Reads from a given file
  ArrayList<String> getReads(String fn) throws Exception
  {
    // ArrayList to store reads
    ArrayList<String> read_list = new ArrayList<String>();

    // Scan file
    Scanner input = new Scanner(new FileInputStream(new File(fn)));

    // While there is another line
    while(input.hasNext())
    {
      String line = input.nextLine();
      // If line is empty or line starts with '>', it is not a read and we continue
      if(line.length() == 0 || line.startsWith(">"))
      {
        continue;
      }
      // Else, add read to the reads ArrayList
      read_list.add(line.toUpperCase());
    }
    // List with all reads
    return read_list;
  }

  // Loads all .fasta/.fa/.fna files in the genomes directory
  List<String> getGenomeFiles(String directory)
  {
    List<String> textFiles = new ArrayList<String>();
    File dir = new File(directory);
    for (File file : dir.listFiles()) {
      if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
        textFiles.add(file.getName());
      }
    }
    return textFiles;
  }

  // Helper function to write results to file
  void saveResults(String filename, String[] readSets, String[] genomeNames, ArrayList<HashSet<Integer>> sketch_hash, int[] totalReads, int[] correctCounts, int[] misCounts, int[] insufCounts, int[] tieCounts) throws Exception
  {
    PrintWriter out = new PrintWriter(new File(filename));
    System.out.println("Saving results...");
    out.println("Genome Sketches:");
    for (int i = 0; i < genomeNames.length; i++)
    {
      out.println(genomeNames[i] + " " + sketch_hash.get(i).size());
    }
    out.println("Read Screening Results:");
    for (int i = 0; i < readSets.length; i++)
    {
      out.println(readSets[i]);
      out.println(totalReads[i] + " " + correctCounts[i] + " " + misCounts[i] + " " + insufCounts[i] + " " + tieCounts[i]);
    }
    out.close();
  }

  // Print the misclassified matrix
  public void printMatrix(int[][] matrix) {
    for (int row = 0; row < matrix.length; row++) {
        for (int col = 0; col < matrix[row].length; col++) {
            System.out.printf("%4d", matrix[row][col]);
        }
        System.out.println();
    }
  }
}
