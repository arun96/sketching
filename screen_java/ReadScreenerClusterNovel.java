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

// Heirachical Clustering Java
import com.apporiented.algorithm.clustering.*;
import com.apporiented.algorithm.clustering.visualization.*;

public class ReadScreenerClusterNovel extends Screener {

  // Main function for screening reads
  // This will take in both a cluster and a screen (the original screen)
  ReadScreenerClusterNovel(Screen sg, ClusterGenerator cg) throws Exception
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

    // Generate full cluster data structure
    System.out.println("Generating Cluster Structure (using downsampling type = " + Settings.DOWNSAMPLE_TYPE +  ", downsampling factor  = " + Settings.DOWNSAMPLE_FACTOR + ")...");
    HashMap<String, HashSet<Integer>> cluster_map = getClusterStructure(cg, sg.sketch_hash);

    // TODO - add option to screen reads independently of the genomes
    // Add ability to track classifications

    System.out.println("Screening Reads - read results will be saved in: " + Settings.READ_LOCATION);

    // TODO - add something to store the classifications

    // Track total reads, assignments, insufficient, and ties
    int[] totalReads = new int[numReadSets];
    int[] insufCounts = new int[numReadSets];
    int[] tieCounts = new int[numReadSets];

    // Loading reads in specified numbers
    if (Settings.IN_CHUNKS){

      // For each organism/element's readset
     for (int r = 0; r < numReadSets; r++){

       boolean fully_read = false;

       int chunk_count = 0;

       int read_start = 0;

       while (!fully_read){

         ArrayList<String> reads = getReadsChunk(readFolder + Settings.READ_SETS[r], chunk_count*Settings.CHUNK, Settings.CHUNK);

         int numReads = reads.size();

         //Screen these reads
         ParallelScreenerClusterNovel ps = new ParallelScreenerClusterNovel(sg, reads, r, read_start, cg, cluster_map);
         ps.run();

         // Counts for this readset
         totalReads[r] += ps.totalReads.intValue();
         insufCounts[r] += ps.insuf.intValue();
         tieCounts[r] += ps.ties.intValue();

         // Print update
         if (Settings.CHUNK_UPDATES) {
           System.out.println(totalReads[r]);
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
       System.out.println(Settings.READ_SETS[r]);
       System.out.println(totalReads[r] + " " + insufCounts[r] + " " + tieCounts[r]);

     }

    // Load file by file
    } else {

      // For each organism/element's readset
     for (int r = 0; r < numGenomes; r++) {

       // Version for automatic loading
       ArrayList<String> reads = getReads(readFolder + Settings.READ_SETS[r]);

       int read_start = 0;

       // For each read in the readset
       int numReads = reads.size();

       ParallelScreenerClusterNovel ps = new ParallelScreenerClusterNovel(sg, reads, r, read_start, cg, cluster_map);
       ps.run();
       // Counts for this readset
       totalReads[r] = ps.totalReads.intValue();
       insufCounts[r] = ps.insuf.intValue();
       tieCounts[r] = ps.ties.intValue();

       // Print summary
       // TODO - fix this for read sets distinct from genomes
       System.out.println(Settings.READ_SETS[r]);
       System.out.println(totalReads[r] + " " + insufCounts[r] + " " + tieCounts[r]);
     }

    }

    // Whole experiment summary
    System.out.println("Printing Results...");
    System.out.println("Number of Reads:");
    System.out.println(Arrays.toString(totalReads));
    System.out.println("Number Not Classified:");
    System.out.println(Arrays.toString(insufCounts));

    //Save results
    saveResults(Settings.OUTPUT_FILE, Settings.READ_SETS, sg.genomeNames, sg.sketch_hash, totalReads, insufCounts, tieCounts);
  }

  // ----- I/O HELPER FUNCTIONS ------

  // Helper function to write results to file
  void saveResults(String filename, String[] readSets, String[] genomeNames, ArrayList<HashSet<Integer>> sketch_hash, int[] totalReads, int[] insufCounts, int[] tieCounts) throws Exception
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
      out.println(totalReads[i] + " " + insufCounts[i] + " " + tieCounts[i]);
    }
    out.close();
  }
}
