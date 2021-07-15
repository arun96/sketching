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


public class ReadClassifierNovel extends Classifier{

  // Read to be processed
  String read;

  // Matches read shares with each sketch
  int[] scores;

  // Data structure to store hashed read
  ArrayList<Integer> read_hashes;

  // Parameters for deciding hashing method
  int window;

  int k;

  // Min number of matches for read to be classified
  int threshold;

  // Read's number and Read Set
  int read_number;
  int read_set;

  // Keep track of the status of the read
  int predicted;
  int score;
  int insufficient;
  int tied;

  boolean filtered_out;

  ReadClassifierNovel(int readSet, String read, int window, int read_number, int sketch_size) {

    scores = new int[sketch_size];

    // this.sketch_hash = sketch_hash;

    this.window = window;

    this.read = read;

    this.k = Settings.K;

    this.threshold = Settings.THRESHOLD;

    this.read_number = read_number;

    this.read_set = readSet;

    this.filtered_out = false;

    // read status
    predicted = 0;
    score = 0;
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

    // Check if tied or insufficient
    if (read_scores[predicted] == 0 || read_scores[predicted] < threshold){
      insufficient++;
    } else {
      if (sameCounts(read_scores, predicted) > 0) {
        tied++;
      }
    }

    // saveReadResults(Settings.READ_LOCATION, read_set, read_number, predicted, score, tied);
    saveReadResults(Settings.READ_LOCATION, read_set, read_number, read_scores);
  }

  // ----- READ LOGGING -----
  // Save read details to a file
  void saveReadResults(String location, int readset, int readnumber, int[] readscores) {
    String filename = location + readset + "_" + readnumber + ".log";
    try {
      PrintWriter out = new PrintWriter(new File(filename));
      out.println(Arrays.toString(readscores));
      out.close();
    } catch (Exception e) {
      e.printStackTrace();
    }
  }

}
