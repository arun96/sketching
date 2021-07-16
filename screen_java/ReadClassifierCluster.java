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


public class ReadClassifierCluster extends Classifier {

  ReadClassifierCluster(String read, int window, int source, int read_number, Cluster cluster, HashMap<String, ArrayList<String>> cluster_sketch_map, HashMap<String, Integer> genome_sketch_map, HashMap<String, Integer> cluster_height_map, int sketch_size){

    scores = new int[sketch_size];

    // this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.source = source;

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
    correct = 0;
    incorrect = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead(Screen sg, HashMap<String, HashSet<Integer>> cluster_map) {

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

    String [] path = screenReadCluster(sg.sketch_hash, cluster_map, cluster, cluster_sketch_map, genome_sketch_map, read_hashes, scores, tie);

    predicted = genome_sketch_map.get(path[path.length - 1]);

    score = (int) scores.get(scores.size() - 1);

    // System.out.println(source + " " + predicted + " " + score);

    // Update counts
    if (predicted == source && score > threshold) {
      // Correctly classified
      correct++;
      if (tie) {
        // Tie was broken correctly
        tied++;
      }
    } else {
      if (score == 0 || score < threshold){
        // Not enough to classify - should not happen often
        insufficient++;
      } else{
        // Misclassified - update counts
        incorrect++;
        if (Settings.TRACK_MISCLASSIFIED){
          System.out.println(source + " " + read_number);
        }
        if (tie){
          // Tie, but this time broken incorrectly
          tied++;
        }
      }
    }

    if (Settings.READ_LOGGING) {
      saveReadResultsCluster(Settings.READ_LOCATION, source, read_number, path, score, predicted, sg.sketch_hash.size());
    }
  }

  // ----- READ LOGGING -----

  // TODO - maybe add more information (number of screens, matches down the path, etc.)
  void saveReadResultsCluster(String location, int source, int readnumber, String[] path, int score, int prediction, int screen_size) {
    String filename = location + source + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(screen_size);
      out.println(Arrays.toString(path));
      out.println(prediction);
      out.println(score);
      out.println(source);
      out.close();
    } catch (Exception e) {
      // Debugging
      // System.out.println(readnumber);
      e.printStackTrace();
    }
  }

}
