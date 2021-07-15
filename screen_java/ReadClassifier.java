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


public class ReadClassifier extends Classifier{

  ReadClassifier(String read, int window, int source, int read_number, int sketch_size){

    scores = new int[sketch_size];

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.source = source;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    this.filtered_out = false;

    // read status
    predicted = 0;
    score = 0;
    correct = 0;
    incorrect = 0;
    insufficient = 0;
    tied = 0;
  }

  // Misclassification matrix can be pieced together from saved read logs
  void classifyRead(ArrayList<HashSet<Integer>> sketch_hash) {

    // System.out.println(source + " " + read_number);

    // TODO - generalize filtering
    // Read Filtering
    if (Settings.FILTER_READS) {
      filtered_out = read_filtering(read, Settings.READ_LENGTH);
      if (filtered_out) {
        // System.out.println("Filtered:" + read_number);
        return;
      }
    }

    // Get read k-mers
    if (window > 0) {
      read_hashes = getAllMinimizers(read, window, k);
    } else {
      read_hashes = getReadKmersHash(read, k);
    }

    int[] read_scores = screenReadHash(sketch_hash, read_hashes);

    predicted = getMaxIndex(read_scores);

    score = read_scores[predicted];
    // System.out.println(score);

    // Update counts
    if (predicted == source && read_scores[predicted] > threshold) {
      // Correctly classified
      correct++;
      if (sameCounts(read_scores, predicted) > 0) {
        // Tie was broken correctly
        tied++;
      }
    } else {
      if (read_scores[predicted] < threshold){
        // Not enough to classify - should not happen often
        insufficient++;
      } else{
        // Misclassified - update counts
        incorrect++;
        if (Settings.TRACK_MISCLASSIFIED){
          System.out.println(source + " " + read_number);
        }
        if (read_scores[predicted] == read_scores[source]){
          // Tie, but this time broken incorrectly
          tied++;
        }
      }
    }

    if (Settings.READ_LOGGING) {
      saveReadResults(Settings.READ_LOCATION, source, read_number, read_scores, source);
    }
  }

  // ----- READ LOGGING -----
  void saveReadResults(String location, int readset, int readnumber, int[] readscores, int s) {
    String filename = location + readset + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(Arrays.toString(readscores));
      out.println(s);
      out.close();
    } catch (Exception e) {

      // Debugging
      // System.out.println(readnumber);

      e.printStackTrace();
    }
  }

}
