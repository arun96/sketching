/*
* Mash Screen
*/

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
import java. util. Collection;

public class ReadScreener {
  // User-specified Variables
  int readLen;
  double readErr;
  int targetMatches;
  String genomeFolder;
  String readFolder;
  String[] genomeNames;
  // String[] genomes;
  String[] readSets;
  int numGenomes;
  ArrayList<ArrayList<String>> sketch;
  int window;

  // Removed the reads list for now - reads are loaded and then replaced, to allow for larger readsets to be tested
  // ArrayList<ArrayList<String>> reads;

  // Variables that are mostly fixed
  //boolean saveFlag;

  public static void main(String[] args) throws Exception
  {

    // UNIVERSAL PARAMETERS

    int k = 21;
    boolean track_minimizers = false;
    boolean track_reads = false;

    // Number of arguments
    int argNum = args.length;

    // Genome Folder - e.g. ./Genomes/
    String gf = args[0];

    // Output directory
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

    // Call the appropriate screening function

    ScreenGenerator sg;

    // Number of target Matches - initialized to zero, then read in
    int tm = 0;

    // TODO - clean this up
    // Uniformly Sampled Screen
    if (argNum == 7 && args[6].equals("u")){
      tm = Integer.parseInt(args[5]);
      UniformScreen us = new UniformScreen(gf, g, rl, re, tm, k, args[6]);
      ReadScreener rs = new ReadScreener(ofn, us, gf, g, r, rf, rl, re, k, track_reads);

    // Minimizer Screen, with calculated window size
    } else if (argNum == 7 && args[6].equals("m")) {
      tm = Integer.parseInt(args[5]);
      MinimizerScreen mzs = new MinimizerScreen(gf, g, rl, re, tm, k, args[6], "c", track_minimizers);
      ReadScreener rs = new ReadScreener(ofn, mzs, gf, g, r, rf, rl, re, k, track_reads);

    // Fixed Size Screen
    } else if (argNum == 7 && args[5].equals("f")) {
      // Screen size
      int ss = Integer.parseInt(args[6]);
      MinHashScreen ms = new MinHashScreen(gf, g, rl, re, k, args[5], ss);
      ReadScreener rs = new ReadScreener(ofn, ms, gf, g, r, rf, rl, re, k, track_reads);

    // Minimizer Screen, with provided window size
    } else if (argNum == 7 && args[5].equals("m")) {
      // Window size
      int ws = Integer.parseInt(args[6]);
      MinimizerScreen mzs = new MinimizerScreen(gf, g, rl, re, k, args[5], ws, "g", track_minimizers);
      ReadScreener rs = new ReadScreener(ofn, mzs, gf, g, r, rf, rl, re, k, track_reads);

    // Default option - MinHash Screen
    } else {
      tm = Integer.parseInt(args[5]);
      MinHashScreen ms = new MinHashScreen(gf, g, rl, re, tm, k);
      ReadScreener rs = new ReadScreener(ofn, ms, gf, g, r, rf, rl, re, k, track_reads);
    }

    System.out.println("Done!");

  }

  // Main function for screening reads
  ReadScreener(String filename, ScreenGenerator sg, String gf, String[] g, String[] r, String rf, int readLength, double readError, int k, boolean track_reads) throws Exception
  {

    boolean READ_TRACKING = track_reads;

    // Store all values
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.readFolder = rf;
    this.genomeNames = g;
    // this.genomes = sg.genomes;
    this.readSets = r;
    this.sketch = sg.sketch;
    this.numGenomes = g.length;
    this.window = sg.window;

    // Print Sketch Sizes
    for (int a = 0; a < numGenomes; a++)
    {
      System.out.println(genomeNames[a] + " " + sketch.get(a).size());
    }

    // For each readset
    System.out.println("Screening Reads...");

    // Initialize count arrays to track accuracy
    int[] totalReads = new int[numGenomes];
    int[] correctCounts = new int[numGenomes];
    int[] misCounts = new int[numGenomes];
    int[] insufCounts = new int[numGenomes];
    int[] tieCounts = new int[numGenomes];

    // Create HashSets of the Sketches
    ArrayList<HashSet<String>> sketch_sets = new ArrayList<HashSet<String>>();
    for (int y = 0; y < numGenomes; y++)
    {
      sketch_sets.add(new HashSet<String>(sketch.get(y)));
    }

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
      //totalReads[rs] = numReads;

      // READ TRACKING
      ArrayList<Integer> match_counts = new ArrayList<Integer>();

      for (int rd = 0; rd < numReads; rd++)
      {
        // Get current read
        String currRead = reads.get(rd);
        if (currRead.contains(";")){
          continue;
        }
        totalReads[rs]++;

        // Get read k-mers
        ArrayList<String> readMers = getReadKmers(currRead);

        // TODO - figure out minimizer read approach
        // ArrayList<String> readMers;
        // if (this.window > 0){
        //   readMers = getReadMinimizers(currRead, window, k);
        // } else {
        //   readMers = getReadKmers(currRead);
        // }

        // Get the number of matches for each read with each genome
        int[] readScores = screenRead(sketch_sets, readMers);

        // Get predicted source
        int pred = getMaxIndex(readScores);

        // READ TRACKING
        match_counts.add(readScores[pred]);

        // Update counts
        if (pred == rs && readScores[pred] > 0) {
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
            // System.out.println(currRead);
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
      // Print update
      System.out.println(readSets[rs]);
      System.out.println(totalReads[rs] + " " + correctCounts[rs] + " " + misCounts[rs] + " " +  insufCounts[rs] + " " + tieCounts[rs]);
      if (READ_TRACKING){
        saveReadMatches(rs, match_counts);
      }
    }

    // Final summary
    System.out.println("Printing Results...");
    System.out.println("Number of Reads:");
    System.out.println(Arrays.toString(totalReads));
    System.out.println("Number Correctly Classified:");
    System.out.println(Arrays.toString(correctCounts));
    System.out.println("Number Incorrectly Classified:");
    System.out.println(Arrays.toString(misCounts));
    System.out.println("Number Not Classified:");
    System.out.println(Arrays.toString(insufCounts));

    printMatrix(miscount_matrix);

    //Save results
    saveResults(filename, readSets, genomeNames, sketch, totalReads, correctCounts, misCounts, insufCounts, tieCounts);
  }


  // ----- HELPER FUNCTIONS FOR SCREENING ------
  // Break Read into Kmers
  ArrayList<String> getReadKmers(String read)
  {
    int readLen = read.length();

    ArrayList<String> readMers = new ArrayList<String>();

    for (int i = 0; i < readLen-21; i++)
    {
      // Slice the kmer, and chose the canonical ones
      String readMer = read.substring(i, i+21);
      String reversedMer = reverseComplement(readMer);
      String selectedMer = "";
      if (readMer.compareTo(reversedMer) > 0){
        selectedMer = reversedMer;
      } else {
        selectedMer = readMer;
      }
      readMers.add(selectedMer);
    }

    return readMers;
  }

  // Screen read kmers against screen
  int[] screenRead(ArrayList<HashSet<String>> sketch_sets, ArrayList<String> readMers){

    // Number of genomes in screen
    int num = sketch.size();

    // Track number of matches
    int[] scores = new int[num];

    // For each sketch in the screen
    for (int i = 0; i < num; i++)
    {
      // Get and store overlap
      int score = getOverlap(sketch_sets.get(i), readMers);
      scores[i] = score;
    }
    return scores;
  }

  // Finds number of shared elements between a sketch and a list of query kmers
  int getOverlap(HashSet<String> sketch_set, ArrayList<String> readMers_list)
  {
    int overlap = 0;

    // DEBUGGING - tracking number of duplicates
    int dupe = 0;

    // DEBUGGING - TRACKING THE NUMBER THAT ARE DUPLICATES
    Set<String> already_counted  = new HashSet<String>();

    for (int i = 0; i < readMers_list.size(); i++) {
      if (sketch_set.contains(readMers_list.get(i))) {
        if (already_counted.contains(readMers_list.get(i))) {
          dupe++;
        } else {
          already_counted.add(readMers_list.get(i));
          overlap++;
        }
      }
    }

    // DEBUGGING - tracking number of duplicates
    //System.out.println(dupe);

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

  // Get the reverse complement of a DNA sequence
  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }


  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Returns a list of minimizers from the reads
  ArrayList<String> getReadMinimizers(String input_read, int window_size, int k) throws Exception
  {
    // Get read lengths
    int read_len = input_read.length();

    // ArrayList of Minimizers
    ArrayList<String> read_minimizers = new ArrayList<String>();

    // Iterate through the read
    int p = 0;
    while (p < read_len - window_size){

      // Get the minimizer
      String current_window = input_read.substring(p, p+window_size);
      String[] min_vals = getMinimizer(current_window, k);
      String current_minimizer = min_vals[0];

      // Get canonical
      String reversed_minimizer = reverseComplement(current_minimizer);
      String selected_minimizer = "";
      if (current_minimizer.compareTo(reversed_minimizer) > 0){
        selected_minimizer = reversed_minimizer;
      } else {
        selected_minimizer = current_minimizer;
      }

      // Add correct minimizer
      read_minimizers.add(selected_minimizer);
    }

    return read_minimizers;

  }

  // Helper function to get the lexicographic minimizer from a given window in a read
  String[] getMinimizer(String window, int k) throws Exception
  {
    String[] ret = new String[2];
    String curr_str = window.substring(0, k);
    String window_min_str = curr_str;
    int min_loc = 0;

    for (int i = 1; i < window.length() - k; i++){

      // Get the next k-1 mer, and the next character
      // curr_str = curr_str.substring(1, k) + window.charAt(i);
      curr_str = window.substring(i, i+k);

      // If it is smaller
      if (window_min_str.compareTo(curr_str) > 0){
        window_min_str = curr_str;
        min_loc = i;
      }
    }
    ret[0] = window_min_str;
    ret[1] = Integer.toString(min_loc);

    return ret;
  }


  // ----- I/O HELPER FUNCTIONS ------
  // Load Reads from file
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

  // Loads all .fasta files in the genomes directory
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
  void saveResults(String filename, String[] readSets, String[] genomeNames, ArrayList<ArrayList<String>> sketch, int[] totalReads, int[] correctCounts, int[] misCounts, int[] insufCounts, int[] tieCounts) throws Exception
  {
    PrintWriter out = new PrintWriter(new File(filename));
    System.out.println("Saving results...");
    out.println("Genome Sketches:");
    for (int i = 0; i < genomeNames.length; i++)
    {
      out.println(genomeNames[i] + " " + sketch.get(i).size());
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


    // ----- READ MATCH TRACKING FUNCTIONS ------

  // Helper function to save number of read matches
  void saveReadMatches(int order, ArrayList<Integer> match_list) throws Exception
  {
    // TODO - re-run this
    String filename = "./tmp/ZYMO/calc/read_matches_hash_double10_" + order + ".txt";

    StringBuilder sb = new StringBuilder();
    for (int i = 0; i <= match_list.size() - 1; i++)
    {
      int num = match_list.get(i);
      sb.append(num);
      sb.append(",");
    }
    String minimizer_string = sb.toString();

    PrintWriter out = new PrintWriter(new File(filename));
    out.println(minimizer_string);
    out.close();
  }



}
