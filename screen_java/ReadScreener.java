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
  int numReadSets;
  ArrayList<HashSet<Integer>> sketch_hash;
  int window;

  // Main function for screening reads
  ReadScreener(ScreenGenerator sg) throws Exception
  {

    // Get the sketch and experiment parameters saved
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.GENOME_FOLDER;
    this.readFolder = Settings.READS_FOLDER;
    this.genomeNames = sg.genomeNames;
    this.readSets = Settings.READ_SETS;
    this.k = Settings.K;

    // The sketch generated by the ScreenGenerator
    this.sketch_hash = sg.sketch_hash;
    this.numGenomes = genomeNames.length;
    this.numReadSets = readSets.length;

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

    // Loading reads in specified numbers
    if (Settings.IN_CHUNKS){

      // For each organism/element's readset
     for (int source = 0; source < numReadSets; source++){

       boolean fully_read = false;

       int chunk_count = 0;

       int read_start = 0;

       while (!fully_read){

         ArrayList<String> reads = getReadsChunk(readFolder + readSets[source], chunk_count*Settings.CHUNK, Settings.CHUNK);

         int numReads = reads.size();

         //Screen these reads
         ParallelScreener ps = new ParallelScreener(sketch_hash, reads, window, source, read_start);
         ps.run();
         // Counts for this readset
         totalReads[source] += ps.totalReads;
         correctCounts[source] += ps.correct.intValue();
         misCounts[source] += ps.mis.intValue();
         insufCounts[source] += ps.insuf.intValue();
         tieCounts[source] += ps.ties.intValue();

         // Print update
         if (Settings.CHUNK_UPDATES) {
           System.out.println(totalReads[source] + " " + correctCounts[source] + " " + misCounts[source] + " " +  insufCounts[source] + " " + tieCounts[source]);
         }

         // We have reached the end of the file
         if (numReads < Settings.CHUNK) {
           fully_read = true;
           break;
         }

         // Update read count
         read_start += numReads;

         // Move to next set of reads
         chunk_count++;

       }

       // Print summary
       System.out.println(readSets[source]);
       System.out.println(totalReads[source] + " " + correctCounts[source] + " " + misCounts[source] + " " +  insufCounts[source] + " " + tieCounts[source]);

     }

    // Load file by file
    } else {

      // For each organism/element's readset
     for (int source = 0; source < numGenomes; source++) {

       // Version for automatic loading
       ArrayList<String> reads = getReads(readFolder + readSets[source]);

       int read_start = 0;

       // For each read in the readset
       int numReads = reads.size();

       ParallelScreener ps = new ParallelScreener(sketch_hash, reads, window, source, read_start);
       ps.run();
       // Counts for this readset
       totalReads[source] = ps.totalReads;
       correctCounts[source] = ps.correct.intValue();
       misCounts[source] = ps.mis.intValue();
       insufCounts[source] = ps.insuf.intValue();
       tieCounts[source] = ps.ties.intValue();

       // Print summary
       System.out.println(readSets[source]);
       System.out.println(totalReads[source] + " " + correctCounts[source] + " " + misCounts[source] + " " +  insufCounts[source] + " " + tieCounts[source]);
     }

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

    // Print misclassification matrix -  this is now done in post processing
    // printMatrix(miscount_matrix);

    //Save results
    saveResults(Settings.OUTPUT_FILE, readSets, genomeNames, sketch_hash, totalReads, correctCounts, misCounts, insufCounts, tieCounts);
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

  // Load a specified section of reads from a given file
  ArrayList<String> getReadsChunk(String fn, int start_point, int num_reads) throws Exception
  {
    // ArrayList to store reads
    ArrayList<String> read_list = new ArrayList<String>();

    // Scan file
    Scanner input = new Scanner(new FileInputStream(new File(fn)));

    // Counters to get the current chunk of reads
    int counter = 0;
    int start = (Settings.READ_LINES*start_point);
    int limit = start + (Settings.READ_LINES*num_reads);

    // While there is another line
    while(input.hasNext())
    {
      String line = input.nextLine();

      // If line has a valid read
      if ((line.length() > 0) && !(line.startsWith(">"))) {
        // If read is in desired range/chunk
        if ((counter >= start) && (counter <= limit)) {
          read_list.add(line.toUpperCase());
        }
      }
      counter++;

      // Break once you have crossed that chunk
      if (counter > limit){
        break;
      }
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
