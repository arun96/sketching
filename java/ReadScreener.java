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

public class ReadScreener {
  // User-specified Variables
  int readLen;
  double readErr;
  int targetMatches;
  String genomeFolder;
  String readFolder;
  String[] genomeNames;
  String[] genomes;
  String[] microbes;
  int numGenomes;
  ArrayList<ArrayList<String>> sketch;

  // Removed the reads list for now - reads are loaded and then replaced, to allow for larger readsets to be tested
  // ArrayList<ArrayList<String>> reads;

  // Variables that are mostly fixed
  //boolean saveFlag;

  public static void main(String[] args) throws Exception
  {

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

    // Target Matches
    int tm = Integer.parseInt(args[5]);

    //Genomes
    /*
    String[] g = {"Bacillus_subtilis_complete_genome.fasta",
        "Cryptococcus_neoformans_draft_genome.fasta",
        "Enterococcus_faecalis_complete_genome.fasta",
        "Escherichia_coli_complete_genome.fasta",
        "Lactobacillus_fermentum_complete_genome.fasta",
        "Listeria_monocytogenes_complete_genome.fasta",
        "Pseudomonas_aeruginosa_complete_genome.fasta",
        "Saccharomyces_cerevisiae_draft_genome.fasta",
        "Salmonella_enterica_complete_genome.fasta",
        "Staphylococcus_aureus_complete_genome.fasta"};
    */

    // List of genomes, when some human genome is included:
    ///*
    String[] g = {"Bacillus_subtilis_complete_genome.fasta",
        "Cryptococcus_neoformans_draft_genome.fasta",
        "Enterococcus_faecalis_complete_genome.fasta",
        "Escherichia_coli_complete_genome.fasta",
        "Lactobacillus_fermentum_complete_genome.fasta",
        "Listeria_monocytogenes_complete_genome.fasta",
        "Pseudomonas_aeruginosa_complete_genome.fasta",
        "Saccharomyces_cerevisiae_draft_genome.fasta",
        "Salmonella_enterica_complete_genome.fasta",
        "Staphylococcus_aureus_complete_genome.fasta",
        "chromosomes/chr18.fa"};
    //*/

    // Microbe Names
    /*
    String[] m = {"Bacillus_subtilis",
        "Cryptococcus_neoformans",
        "Enterococcus_faecalis",
        "Escherichia_coli",
        "Lactobacillus_fermentum",
        "Listeria_monocytogenes",
        "Pseudomonas_aeruginosa",
        "Saccharomyces_cerevisiae",
        "Salmonella_enterica",
        "Staphylococcus_aureus"};
    */

    // Microbe Names + Human
    ///*
    String[] m = {"Bacillus_subtilis",
        "Cryptococcus_neoformans",
        "Enterococcus_faecalis",
        "Escherichia_coli",
        "Lactobacillus_fermentum",
        "Listeria_monocytogenes",
        "Pseudomonas_aeruginosa",
        "Saccharomyces_cerevisiae",
        "Salmonella_enterica",
        "Staphylococcus_aureus",
        "Human/chr18"};
    //*/

    // Call the get screen function
    System.out.println("Generating Screen...");
    ScreenGenerator sg = new ScreenGenerator(gf, g, rl, re, tm, 21);

    ReadScreener rs = new ReadScreener(ofn, sg, gf, g, m, rf, rl, re);

    //rs.saveResults();
  }

  ReadScreener(String filename, ScreenGenerator sg, String gf, String[] g, String[] m, String rf, int readLength, double readError) throws Exception
  {
    this.readLen = readLength;
    this.readErr = readError;
    this.genomeFolder = gf;
    this.readFolder = rf;
    this.genomeNames = g;
    this.genomes = sg.genomes;
    this.microbes = m;
    this.sketch = sg.sketch;
    this.numGenomes = g.length;

    // Sketch Sizes
    for (int a = 0; a < numGenomes; a++)
    {
      System.out.println(sketch.get(a).size());
    }
    // System.out.println(sketch.get(0));

    // For each readset
    System.out.println("Screening Reads...");

    // Initialize count arrays
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

    for (int rs = 0; rs < numGenomes; rs++)
    {
      // Initialize counts for this readset
      correctCounts[rs] = 0;
      misCounts[rs] = 0;
      insufCounts[rs] = 0;
      totalReads[rs] = 0;
      tieCounts[rs] = 0;

      // Load the readset - not saved any more
      ArrayList<String> reads = getReads(readFolder + microbes[rs] + ".fasta");

      // For each read in the readset
      int numReads = reads.size();
      //totalReads[rs] = numReads;

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

        // Get the number of matches for each read with each genome
        int[] readScores = screenRead(sketch_sets, readMers);

        // Get predicted source
        int pred = getMaxIndex(readScores);

        // Update counts
        if (pred == rs && readScores[pred] > 0) {
          correctCounts[rs]++;
          if (sameCounts(readScores, pred) > 0) {
            tieCounts[rs]++;
          }
        } else {
          if (readScores[pred] == 0){
            insufCounts[rs]++;
            // System.out.println(currRead);
          } else{
            misCounts[rs]++;
            miscount_matrix[rs][pred]++;
            if (readScores[pred] == readScores[rs]){
              tieCounts[rs]++;
            }
          }
        }
      }
      System.out.println(microbes[rs]);
      System.out.println(totalReads[rs] + " " + correctCounts[rs] + " " + misCounts[rs] + " " +  insufCounts[rs] + " " + tieCounts[rs]);
    }

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
    saveResults(filename, microbes, sketch, totalReads, correctCounts, misCounts, insufCounts, tieCounts);

    System.out.println("Done!");
  }

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

    return read_list;
  }

  ArrayList<String> getReadKmers(String read)
  {
    int readLen = read.length();

    ArrayList<String> readMers = new ArrayList<String>();

    for (int i = 0; i < readLen-21; i++)
    {
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

  int[] screenRead(ArrayList<HashSet<String>> sketch_sets, ArrayList<String> readMers){

    int num = sketch.size();
    int[] scores = new int[num];
    for (int i = 0; i < num; i++)
    {
      int score = getOverlap(sketch_sets.get(i), readMers);
      scores[i] = score;
    }
    return scores;
  }

  int getOverlap(HashSet<String> sketch_set, ArrayList<String> readMers_list)
  {
    // New Set Approach
    int overlap = 0;
    for (int i = 0; i < readMers_list.size(); i++) {
      if (sketch_set.contains(readMers_list.get(i))) {
        overlap++;
      }
    }
    return overlap;
  }

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

  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }

  void saveResults(String filename, String[] microbes, ArrayList<ArrayList<String>> sketch, int[] totalReads, int[] correctCounts, int[] misCounts, int[] insufCounts, int[] tieCounts) throws Exception
  {
    PrintWriter out = new PrintWriter(new File(filename));
    System.out.println("Saving results...");
    for (int i = 0; i < microbes.length; i++)
    {
      out.println(microbes[i] + " " + sketch.get(i).size());
      out.println(totalReads[i] + " " + correctCounts[i] + " " + misCounts[i] + " " + insufCounts[i] + " " + tieCounts[i]);
    }
    out.close();
  }

  public void printMatrix(int[][] matrix) {
    for (int row = 0; row < matrix.length; row++) {
        for (int col = 0; col < matrix[row].length; col++) {
            System.out.printf("%4d", matrix[row][col]);
        }
        System.out.println();
    }
  }

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

}
