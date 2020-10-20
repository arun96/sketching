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
import java.util.Collection;

public class MinimizerScreen {
  // User-specified Variables
  int targetMatches;
  int readLen;
  double readErr;
  String genomeFolder;
  int k;
  int numGenomes;
  int[] sketch_size;
  String[] genomeNames;
  ArrayList<ArrayList<String>> sketch;
  String screenType;
  int window;

  public static void main(String[] args) throws Exception
  {

    // Constructor

  }

  // Screen Generator for minimzer-based approach, with calculated windowsize
  MinimizerScreen(String gf, String[] g, int readLength, double readError, int tm, int kmer, String mini, String calc, boolean track) throws Exception
  {
    // TRUE if we want to save minimizer positions
    boolean MINIMIZER_TRACKING = track;

    // Adjust
    double multiplier = 2.0;


    System.out.println("Generating Minimizer-Based Screen with calculated Window Size...");
    this.screenType = "v";

    // Store variables
    this.targetMatches = tm;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;

    // Get genome lengths
    int[] genomeLengths = new int[numGenomes];

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + g[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
      // System.out.println(genomeLengths[i]);
    }

    // Compute sketch size
    sketch_size = new int[numGenomes];
    int[] window_sizes = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
      window_sizes[j] = (int) ((genomeLengths[j]/sketch_size[j])*multiplier);
    }

    // set the window size - agnostic of genomes
    this.window = window_sizes[0];

    // Get sketch using genomes and window sizes

    sketch = new ArrayList<ArrayList<String>>();

    for (int x = 0; x < numGenomes; x++)
    {
      // System.out.println(x);

      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // MINIMIZER TRACKING - arraylist for minimizer positions
      ArrayList<Integer> minimizer_positions = new ArrayList<Integer>();
      int prev_mininimizer_pos = 0;
      int current_minimizer_pos = 0;

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();

      // Get the window size for this sketch
      int window_size = window_sizes[x];

      int p = 0;

      while (p < genomeLengths[x] - window_size) {

        // Get the current window
        String current_window = genomes[x].substring(p, p+window_size);

        // Get the minimzer for this window
        // LEXICOGRAPHIC
        // String[] min_vals = getMinimizer(current_window, k);
        // HASH BASED
        String[] min_vals = getMinimizerHash(current_window, k);
        String current_minimizer = min_vals[0];

        //MINIMIZER TRACKING - Tracking the minimizer position
        if (MINIMIZER_TRACKING) {

          // Get the minimizer position
          current_minimizer_pos = (p + Integer.parseInt(min_vals[1]));
          // If it's a new minimizer
          if (current_minimizer_pos != prev_mininimizer_pos){
            // Update the prev
            prev_mininimizer_pos = current_minimizer_pos;
            // Save the position
            minimizer_positions.add(current_minimizer_pos);

          // Same minimizer as before
          } else {
            prev_mininimizer_pos = prev_mininimizer_pos;
          }
        }

        // Get canonically smaller of minimzer and its reverse complement
        String reversed_minimizer = reverseComplement(current_minimizer);
        String selected_minimizer = "";


        // HASH CODE BASED APPROACH - NOT CURRENTLY USED
        // if (current_minimizer.hashCode() > reversed_minimizer.hashCode()){
        //   selected_minimizer = reversed_minimizer;
        // } else {
        //   selected_minimizer = current_minimizer;
        // }


        if (current_minimizer.compareTo(reversed_minimizer) > 0){
          selected_minimizer = reversed_minimizer;
        } else {
          selected_minimizer = current_minimizer;
        }

        // If we already have this minimizer
        if (sketch_set.contains(selected_minimizer)){
          p++;
        }
        // If not, add it to the sketch
        else {
          sketch_set.add(selected_minimizer);
          sketch.get(x).add(selected_minimizer);
          p++;
        }
      }

      // MINIMIZER TRACKING - save the minimizer positions to a file
      if (MINIMIZER_TRACKING){
        saveMinimizers(x, minimizer_positions);
      }

    }
  }

  // Screen Generator for minimzer-based approach, with specified window size
  MinimizerScreen(String gf, String[] g, int readLength, double readError, int kmer, String mini, int windowSize, String given, boolean track) throws Exception
  {
    // TRUE if we want to save minimizer positions
    boolean MINIMIZER_TRACKING = track;

    System.out.println("Generating Minimizer-Based Screen with specified Window Size...");

    System.out.println(windowSize);

    // Store variables
    this.targetMatches = 0;
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.k = kmer;
    this.numGenomes = g.length;
    this.genomeNames = g;
    this.window = windowSize;

    // Get genome lengths
    int[] genomeLengths = new int[numGenomes];

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + g[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
      // System.out.println(genomeLengths[i]);
    }

    // Compute sketch size
    int[] window_sizes = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      window_sizes[j] = windowSize;
    }

    // Get sketch using genomes and window sizes

    sketch = new ArrayList<ArrayList<String>>();

    for (int x = 0; x < numGenomes; x++)
    {

      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // MINIMIZER TRACKING - arraylist for minimizer positions
      ArrayList<Integer> minimizer_positions = new ArrayList<Integer>();
      int prev_mininimizer_pos = 0;
      int current_minimizer_pos = 0;

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();

      // Get the window size for this sketch
      int window_size = window_sizes[x];

      int p = 0;
      while (p < genomeLengths[x] - window_size) {

        // Get the current window
        String current_window = genomes[x].substring(p, p+window_size);

        // Get the minimzer for this window

        // LEXICOGRAPHIC
        // String[] min_vals = getMinimizer(current_window, k);

        // HASH BASED
        String[] min_vals = getMinimizerHash(current_window, k);
        String current_minimizer = min_vals[0];

        //MINIMIZER TRACKING - Tracking the minimizer position
        if (MINIMIZER_TRACKING) {

          // Get the minimizer position
          current_minimizer_pos = (p + Integer.parseInt(min_vals[1]));
          // If it's a new minimizer
          if (current_minimizer_pos != prev_mininimizer_pos){
            // Update the prev
            prev_mininimizer_pos = current_minimizer_pos;
            // Save the position
            minimizer_positions.add(current_minimizer_pos);

          // Same minimizer as before
          } else {
            prev_mininimizer_pos = prev_mininimizer_pos;
          }
        }

        // Get canonically smaller of minimzer and its reverse complement
        String reversed_minimizer = reverseComplement(current_minimizer);
        String selected_minimizer = "";
        if (current_minimizer.compareTo(reversed_minimizer) > 0){
          selected_minimizer = reversed_minimizer;
        } else {
          selected_minimizer = current_minimizer;
        }

        // If we already have this minimizer
        if (sketch_set.contains(selected_minimizer)){
          p++;
        }
        // If not, add it to the sketch
        else {
          sketch_set.add(selected_minimizer);
          sketch.get(x).add(selected_minimizer);
          p++;
        }
      }

      // MINIMIZER TRACKING - save the minimizer positions to a file
      if (MINIMIZER_TRACKING){
        saveMinimizers(x, minimizer_positions);
      }
    }
  }


  // ----- UTILITY FUNCTIONS ------
  // Compute the sketch size, given a genome, read lengths, error rates and # of target matches
  int getSketchSize(int genomeLength, int readLength, double readError, int targetMatches, int k){
    double unaffectedChance = 1.0 - readError;
    // System.out.println(unaffectedChance);
    double multiplier = Math.pow(unaffectedChance, k);
    // System.out.println(multiplier);

    double genomeLengthD = (double) genomeLength;
    double readLengthD = (double) readLength;
    double targetMatchesD = (double) targetMatches;

    double numerator = (targetMatchesD * genomeLengthD);
    // System.out.println(numerator);
    double denominator = (readLengthD * multiplier);
    // System.out.println(denominator);

    double num_sketches = numerator/denominator;

    return (int) num_sketches;
  }

  // Generates the reverse complement of a DNA sequence
  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  // ----- I/O HELPER FUNCTIONS ------
  // Load genome from a given file
  String getGenome(String fn) throws Exception
  {
    Scanner input = new Scanner(new FileInputStream(new File(fn)));
    StringBuilder sb = new StringBuilder("");
    while(input.hasNext())
    {
      String line = input.nextLine();
      if(line.length() == 0 || line.startsWith(">"))
      {
        continue;
      }
      sb.append(line.toUpperCase());
    }

    String s = sb.toString();

    return s;
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


  // ----- MINIMIZER HELPER FUNCTIONS ------
  // Helper function to get the lexicographic minimizer from a given window
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

  // Helper function to get the hash based minimizer from a given window
  String[] getMinimizerHash(String window, int k) throws Exception
  {
    String[] ret = new String[2];
    String curr_str = window.substring(0, k);
    String window_min_str = curr_str;
    int min_loc = 0;
    int curr_min = curr_str.hashCode();
    int window_min = curr_min;

    for (int i = 1; i < window.length() - k; i++){

      // Get the next k-1 mer, and the next character
      // curr_str = curr_str.substring(1, k) + window.charAt(i);
      curr_str = window.substring(i, i+k);
      curr_min = curr_str.hashCode();

      // If it is smaller
      if (curr_min < window_min) {
        // Update
        window_min_str = curr_str;
        window_min = curr_min;
        min_loc = i;
      }
    }
    ret[0] = window_min_str;
    ret[1] = Integer.toString(min_loc);

    return ret;
  }


  // ----- MINIMIZER TRACKING HELPER FUNCTIONS ------
  // Helper Function to save the minimizers from the genome
  void saveMinimizers(int order, ArrayList<Integer> minimizer_list) throws Exception
  {
    String filename = "./tmp/ZYMO/calc/minimizers_hash_double10_" + order + ".txt";

    // String minimizer_string = String.join(",", minimizer_list);
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i <= minimizer_list.size() - 1; i++)
    {
      int num = minimizer_list.get(i);
      sb.append(num);
      sb.append(",");
    }
    String minimizer_string = sb.toString();

    PrintWriter out = new PrintWriter(new File(filename));
    out.println(minimizer_string);
    out.close();
  }
}
