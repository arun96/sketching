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
  int k;
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
    String genome_folder = args[0];

    // Output file
    String ofn = args[1];

    // Read Folder
    String read_folder = args[2];

    // Read Length
    int read_length = Integer.parseInt(args[3]);

    // Read error
    double read_error = Double.parseDouble(args[4]);

    // Get a sorted list list of all genomes in the input folder
    List<String> genomes_list = new ArrayList<String>();
    File dir = new File(genome_folder);
    for (File file : dir.listFiles()) {
      if (file.getName().endsWith((".fasta")) || file.getName().endsWith((".fna")) || file.getName().endsWith((".fa"))) {
        genomes_list.add(file.getName());
      }
    }
    Collections.sort(genomes_list);
    String[] genomes = new String[genomes_list.size()];
    genomes = genomes_list.toArray(genomes);

    // Get a sorted list of all readsets in the input folder
    List<String> reads_list = new ArrayList<String>();
    File readDir = new File(read_folder);
    for (File readFile : readDir.listFiles()) {
      if (readFile.getName().endsWith((".fasta")) || readFile.getName().endsWith((".fna")) || readFile.getName().endsWith((".fa"))) {
        reads_list.add(readFile.getName());
      }
    }
    Collections.sort(reads_list);
    String[] read_sets = new String[reads_list.size()];
    read_sets = reads_list.toArray(read_sets);

    // ----- SCREENING AND CLASSIFICATION -----
    // Number of target Matches - initialized to zero, then read in
    int tm = 0;

    // Hash options
    // TODO - Add more, and update "getHashName" in ScreenGenerator
    List<String> hashTypes = Arrays.asList("h", "mmh3", "mmh3_128");
    // Default is built-in hashcode
    String hashType = "h";

    // Absolute base option - MinHash Screen
    if (argNum == 6){
      tm = Integer.parseInt(args[5]);
      MinHashScreen ms = new MinHashScreen(genome_folder, genomes, read_length, read_error, tm, k, hashType);
      ReadScreener rs = new ReadScreener(ofn, ms, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

    // MinHash Screen with specified hash function
    } else if ((argNum == 7) && hashTypes.contains(args[6])){
      tm = Integer.parseInt(args[5]);
      hashType = args[6];
      MinHashScreen ms = new MinHashScreen(genome_folder, genomes, read_length, read_error, tm, k, hashType);
      ReadScreener rs = new ReadScreener(ofn, ms, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

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
        UniformScreen us = new UniformScreen(genome_folder, genomes, read_length, read_error, tm, k, args[6], hashType);
        ReadScreener rs = new ReadScreener(ofn, us, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

      // Minimizer Screen, with calculated window size
      } else if (args[6].equals("m")) {
        tm = Integer.parseInt(args[5]);
        MinimizerScreen mzs = new MinimizerScreen(genome_folder, genomes, read_length, read_error, tm, k, args[6], hashType);
        ReadScreener rs = new ReadScreener(ofn, mzs, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

      // Minimizer Screen, with provided window size
      } else if (args[5].equals("m")) {
        // Window size
        int ws = Integer.parseInt(args[6]);
        MinimizerScreen mzs = new MinimizerScreen(genome_folder, genomes, read_length, read_error, k, args[5], ws, hashType);
        ReadScreener rs = new ReadScreener(ofn, mzs, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

      // Fixed Size MinHash Screen
      } else if (args[5].equals("f")) {
        // Screen size
        int ss = Integer.parseInt(args[6]);
        MinHashScreen ms = new MinHashScreen(genome_folder, genomes, read_length, read_error, k, args[5], ss, hashType);
        ReadScreener rs = new ReadScreener(ofn, ms, genome_folder, genomes, read_sets, read_folder, read_length, read_error, k, track_reads, hashType);

      // Invalid parameters
      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }

    }
    System.out.println("Done!");
  }

  // Main function for screening reads
  ReadScreener(String filename, ScreenGenerator sg, String genome_folder, String[] g, String[] r, String read_folder, int read_length, double read_error, int k, boolean track_reads, String hashType) throws Exception
  {
    // Minimum number of matches to classify a read
    int threshold = 2;

    // Number of threads to use - TODO make this a parameter
    int num_threads = 4;

    // Get the sketch and experiment parameters saved
    this.readLen = read_length;
    this.readErr = read_error;
    this.genomeFolder = genome_folder;
    this.readFolder = read_folder;
    this.genomeNames = g;
    this.readSets = r;
    this.k = k;
    // The sketch generated by the ScreenGenerator
    this.sketch_hash = sg.sketch_hash;
    this.numGenomes = g.length;
    // 0 if MinHash or uniform, >0 if Minimizer-based
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

      // Version for automatic loading
      ArrayList<String> reads = getReads(readFolder + readSets[rs]);

      // For each read in the readset
      int numReads = reads.size();

      //TODO - move to parallel
      ParallelScreener ps = new ParallelScreener(sketch_hash, reads, hashType, k, num_threads, window, rs, threshold);
      ps.run();
      // Counts for this readset
      totalReads[rs] = ps.totalReads;
      correctCounts[rs] = ps.correct.intValue();
      misCounts[rs] = ps.mis.intValue();
      insufCounts[rs] = ps.insuf.intValue();
      tieCounts[rs] = ps.ties.intValue();

      // for (int rd = 0; rd < numReads; rd++)
      // {
      //   // Get current read
      //   String currRead = reads.get(rd);
      //   totalReads[rs]++;
      //
      //   ReadClassifier rc = new ReadClassifier(sketch_hash, currRead, hashType, k, window, threshold, rs);
      //   rc.classifyRead();
      //
      //   correctCounts[rs] = correctCounts[rs] + rc.correct;
      //   misCounts[rs] = misCounts[rs] + rc.incorrect;
      //   tieCounts[rs] = tieCounts[rs] + rc.tied;
      //   insufCounts[rs] = insufCounts[rs] + rc.insufficient;
      // }
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

    // Print misclassification matrix - TODO, fix this.
    // printMatrix(miscount_matrix);

    //Save results
    saveResults(filename, readSets, genomeNames, sketch_hash, totalReads, correctCounts, misCounts, insufCounts, tieCounts);
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
