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

public class MinimizerScreen extends ScreenGenerator{

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

}
