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

public class ScreenGenerator {
  // User-specified Variables
  int targetMatches;
  int readLen;
  double readErr;
  String genomeFolder;
  int k;
  int numGenomes;
  int[] sketch_size;
  String[] genomeNames;
  String[] genomes;
  ArrayList<ArrayList<String>> sketch;

  public static void main(String[] args) throws Exception
  {

    // Genome Folder - e.g. ./Genomes/
    String gf = args[0];

    // List of genomes to use - maybe parameterize this later?

    /*
    String[] genomeNames = {"Bacillus_subtilis_complete_genome.fasta","Cryptococcus_neoformans_draft_genome.fasta",
        "Enterococcus_faecalis_complete_genome.fasta",
        "Escherichia_coli_complete_genome.fasta",
        "Lactobacillus_fermentum_complete_genome.fasta",
        "Listeria_monocytogenes_complete_genome.fasta",
        "Pseudomonas_aeruginosa_complete_genome.fasta",
        "Saccharomyces_cerevisiae_draft_genome.fasta",
        "Salmonella_enterica_complete_genome.fasta",
        "Staphylococcus_aureus_complete_genome.fasta"};
    */

    // List of genomes, when the human genome is included:
    ///*
    String[] genomeNames = {"Bacillus_subtilis_complete_genome.fasta",
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

    // Read Length
    int rl = Integer.parseInt(args[2]);

    // Read error
    double re = Double.parseDouble(args[3]);

    // Target Matches
    int tm = Integer.parseInt(args[4]);

    // Generate the screen
    ScreenGenerator sg = new ScreenGenerator(gf, genomeNames, rl, re, tm, 21);

  }

  ScreenGenerator(String gf, String[] g, int readLength, double readError, int tm, int kmer) throws Exception
  {
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
    this.genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + g[i]);
      genomes[i] = gnm;
      genomeLengths[i] = gnm.length();
      //System.out.println(genomeLengths[i]);
    }

    // Compute sketch size
    sketch_size = new int[numGenomes];
    for (int j = 0; j < numGenomes; j++)
    {
      sketch_size[j] = getSketchSize(genomeLengths[j], readLen, readErr, targetMatches, k);
    }

    // Get sketch using genomes and sketch sizes
    // TODO - Turn this into a hash function. For now, it's just random start positions.
    Random r = new Random();

    sketch = new ArrayList<ArrayList<String>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Row corresponding to this genome
      sketch.add(new ArrayList<String>());

      // Set for tracking which kmers are in sketch
      Set<String> sketch_set  = new HashSet<String>();

      // For each sketch in a given genome (there are a total of sketch_size[x] of them)
      int p = 0;
      while (p < sketch_size[x]){
        // Get k-mer
        int start = r.nextInt(genomeLengths[x]-k+1);
        int end = start + k;
        String mer = genomes[x].substring(start, end);

        // Get canonical
        String reversed_mer = reverseComplement(mer);
        String selected_mer = "";
        if (mer.compareTo(reversed_mer) > 0){
          selected_mer = reversed_mer;
        } else {
          selected_mer = mer;
        }

        // If it is already in the sketch
        if (sketch_set.contains(selected_mer)){
          continue;
        // Add to sketch, increment count
        } else {
          sketch_set.add(selected_mer);
          sketch.get(x).add(selected_mer);
          p++;
        }
      }
    }
  }

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

  String reverseComplement(String sequence)
  {
    String reversed_tmp = sequence.replace("A", "t").replace("T", "a").replace("C", "g").replace("G", "c").toUpperCase();
    String reversed = new StringBuffer(reversed_tmp).reverse().toString();
    return reversed;
  }
}
