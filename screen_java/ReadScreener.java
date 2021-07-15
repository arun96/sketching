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

public class ReadScreener extends Screener {

  // Main function for screening reads
  ReadScreener(Screen sg) throws Exception
  {

    // Get the sketch and experiment parameters saved
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.GENOME_FOLDER;
    this.readFolder = Settings.READS_FOLDER;
    this.k = Settings.K;
    this.numGenomes = sg.genomeNames.length;
    this.numReadSets = Settings.READ_SETS.length;

    // 0 if MinHash or uniform, >0 if Minimizer-based
    this.window = sg.window;


    // Print Sketch Sizes
    for (int a = 0; a < numGenomes; a++)
    {
      System.out.println(sg.genomeNames[a] + " " + sg.sketch_hash.get(a).size());
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

         ArrayList<String> reads = getReadsChunk(readFolder + Settings.READ_SETS[source], chunk_count*Settings.CHUNK, Settings.CHUNK);

         int numReads = reads.size();

         //Screen these reads
         ParallelScreener ps = new ParallelScreener(sg, reads, source, read_start);
         ps.run();
         // Counts for this readset
         totalReads[source] += ps.totalReads.intValue();
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
       System.out.println(Settings.READ_SETS[source]);
       System.out.println(totalReads[source] + " " + correctCounts[source] + " " + misCounts[source] + " " +  insufCounts[source] + " " + tieCounts[source]);

     }

    // Load file by file
    } else {

      // For each organism/element's readset
     for (int source = 0; source < numGenomes; source++) {

       // Version for automatic loading
       ArrayList<String> reads = getReads(readFolder + Settings.READ_SETS[source]);

       int read_start = 0;

       // For each read in the readset
       int numReads = reads.size();

       ParallelScreener ps = new ParallelScreener(sg, reads, source, read_start);
       ps.run();
       // Counts for this readset
       totalReads[source] = ps.totalReads.intValue();
       correctCounts[source] = ps.correct.intValue();
       misCounts[source] = ps.mis.intValue();
       insufCounts[source] = ps.insuf.intValue();
       tieCounts[source] = ps.ties.intValue();

       // Print summary
       System.out.println(Settings.READ_SETS[source]);
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

    //Save results
    saveResults(Settings.OUTPUT_FILE, Settings.READ_SETS, sg.genomeNames, sg.sketch_hash, totalReads, correctCounts, misCounts, insufCounts, tieCounts);
  }

  // ----- I/O HELPER FUNCTIONS ------

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
}
