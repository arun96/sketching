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
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

// Heirachical Clustering Java
import com.apporiented.algorithm.clustering.*;
import com.apporiented.algorithm.clustering.visualization.*;


public class ReadClassifierClusterNovel extends Classifier {


  ReadClassifierClusterNovel(String read, int window, int readSet, int read_number, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, HashMap<String, Integer> cluster_height_map, int sketch_size){

    scores = new int[sketch_size];

    // this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.read_set = readSet;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    this.filtered_out = false;

    // Cluster Information
    this.cluster = cluster;
    this.cluster_sketch_map = cluster_sketch_map;
    this.genome_sketch_map = genome_sketch_map;
    this.cluster_height_map = cluster_height_map;

    // read status
    predicted = 0;
    score = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead(ArrayList<HashSet<Integer>> sketch_hash, HashMap<String, HashSet<Integer>> cluster_map) {

    // TODO - generalize filtering
    // Read Filtering
    if (Settings.FILTER_READS) {
      filtered_out = read_filtering(read, Settings.READ_LENGTH);
      if (filtered_out) {
        return;
      }
    }

    // Get read k-mers
    if (window > 0) {
      read_hashes = getAllMinimizers(read, window, k);
    } else {
      read_hashes = getReadKmersHash(read, k);
    }

    boolean tie = false;
    ArrayList<Integer> scores = new ArrayList<Integer>();

    String [] path = screenReadCluster(sketch_hash, cluster_map, cluster, cluster_sketch_map, genome_sketch_map, read_hashes, scores, tie);

    predicted = genome_sketch_map.get(path[path.length - 1]);

    score = (int) scores.get(scores.size() - 1);

    // System.out.println(predicted + " " + score + " " + read_set + " " + read_number);

    // Check if insufficient or tied
    if (score == 0 || score < threshold) {
      // Not enough to classify - should not happen often
      insufficient++;
    }
    if (tie) {
      tied++;
    }

    if (Settings.READ_LOGGING) {
      saveReadResultsCluster(Settings.READ_LOCATION, read_set, read_number, path, score, predicted, sketch_hash.size());
    }
  }

  // ----- HELPER FUNCTIONS FOR SCREENING ------

  // Screen read kmers against clustered approach
  String[] screenReadCluster(ArrayList<HashSet<Integer>> sketch, HashMap<String, HashSet<Integer>> cluster_map, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, ArrayList<Integer> readMers, ArrayList<Integer> scores, boolean tie){

    int pred = 0;
    int selected = 0;
    int latest_score = 0;
    boolean leaf = false;
    ArrayList<String> path_list = new ArrayList<String>();

    while(true) {

      if (cluster.isLeaf()) {
        // path_list.add(cluster.getName());
        break;
      }

      List<Cluster> children = cluster.getChildren();


      int[] matches = new int[children.size()];

      // For each child
      for (int i = 0; i < children.size(); i++) {

        String n = children.get(i).getName();
        // System.out.println("Current Cluster:");
        // System.out.println(n);

        if (children.get(i).isLeaf()) {

          int si = genome_sketch_map.get(n);
          int s = getOverlapHash(sketch.get(si), readMers);

          matches[i] = s;

        } else {

          int ss = getOverlapHash(cluster_map.get(n), readMers);
          matches[i] = ss;
        }
      }

      // Get the cluster to descend into
      selected = getMaxIndex(matches);

      // Checking for tie
      if (sameCounts(matches, selected) > 0) {
        tie = true;
      }

      // Descend into cluster and update
      cluster = children.get(selected);
      path_list.add(cluster.getName());
      scores.add(matches[selected]);
    }

    // Return the path that was taken
    //int[] path_list = path.stream().mapToInt(i -> i).toArray();
    String[] path = path_list.toArray(new String[0]);
    return path;
  }

  // ----- READ LOGGING -----

  // TODO - finalize
  void saveReadResultsCluster(String location, int read_set, int readnumber, String[] path, int score, int prediction, int screen_size) {
    String filename = location + read_set + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(screen_size);
      out.println(Arrays.toString(path));
      out.println(prediction);
      out.println(score);
      out.close();
    } catch (Exception e) {
      // Debugging
      // System.out.println(readnumber);
      e.printStackTrace();
    }
  }

}
