/*
* Generate the screen for input genomes
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
import java.util.Collection;

public class ScreenGenerator {
  // User-specified Variables
  public int targetMatches;
  public int readLen;
  public double readErr;
  public String genomeFolder;
  public int k;
  public int numGenomes;
  public int[] sketch_size;
  public String[] genomeNames;
  public ArrayList<ArrayList<String>> sketch;
  public String screenType;
  public int window;

  public static void main(String[] args) throws Exception
  {

    // Constructor

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
