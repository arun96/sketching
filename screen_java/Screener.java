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

public class Screener {

  public int k;
  public int readLen;
  public double readErr;
  public int targetMatches;
  public String genomeFolder;
  public String readFolder;
  public String[] genomeNames;
  public String[] readSets;
  public int numGenomes;
  public int numReadSets;
  public int window;

  public static void main(String[] args) throws Exception
  {
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

  // // Print the misclassified matrix
  // public void printMatrix(int[][] matrix) {
  //   for (int row = 0; row < matrix.length; row++) {
  //       for (int col = 0; col < matrix[row].length; col++) {
  //           System.out.printf("%4d", matrix[row][col]);
  //       }
  //       System.out.println();
  //   }
  // }

  // ------ CLUSTERING -------

  HashMap<String, HashSet<Integer>> getClusterStructure(ClusterGenerator cg, ArrayList<HashSet<Integer>> sketch) {

    HashMap<String, HashSet<Integer>> map = new HashMap<String, HashSet<Integer>>();

    // Iterate through all clusters
    for ( String key : cg.cluster_sketch_map.keySet() ) {

      // Get list of all sketches under this cluster
      ArrayList<String> sketch_list = cg.cluster_sketch_map.get(key);

      // Get the indices of all sketches
      ArrayList<Integer> sketch_indices_list = new ArrayList<Integer>();
      for (int j = 0; j < sketch_list.size(); j++) {
        String sn = sketch_list.get(j);
        sketch_indices_list.add(cg.genome_sketch_map.get(sn));
      }

      // Get the heights
      int height = cg.cluster_height_map.get(key);

      // Get the downsampled sketches, combined for this cluster
      HashSet<Integer> cluster_sketches = new HashSet<Integer>();
      for (int k = 0; k < sketch_indices_list.size(); k++) {
        // For each sketch, get the whole sketch
        List<Integer> sketch_vals = new LinkedList<Integer>(sketch.get(sketch_indices_list.get(k)));

        // Downsampling
        int downsample_factor = 1;

        // No downsampling
        if (Settings.DOWNSAMPLE_TYPE.equals("n")) {
          downsample_factor = 1;
        // Constant
        } else if (Settings.DOWNSAMPLE_TYPE.equals("c")) {
          downsample_factor = Settings.DOWNSAMPLE_FACTOR;
        // Height based
        } else if (Settings.DOWNSAMPLE_TYPE.equals("h")) {
          downsample_factor = (int) (Math.pow(Settings.DOWNSAMPLE_FACTOR,height));
        }

        // Compute size/sketch based on that
        int downsample_size = (int) (sketch_vals.size() / downsample_factor);

        // Downsample the sketch
        Collections.shuffle(sketch_vals);
        Set<Integer> downsampled_sketch = new HashSet<Integer>(sketch_vals.subList(0, downsample_size));

        // Add it
        cluster_sketches.addAll(downsampled_sketch);
      }

      // Add to the map
      map.put(key, cluster_sketches);
    }
    return map;
  }
}
