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
import java.util.Random;

// GUAVA
import com.google.common.hash.*;
import com.google.common.hash.Hashing;
import com.google.common.hash.HashFunction;

// CLUSTERING
import com.clust4j.algo.*;


public class ClusterScreen{

  ArrayList<HashSet<Integer>> cluster_screen;

  ClusterScreen(ScreenGenerator screen, int[] assignments) throws Exception{

    System.out.println("Clustering generated sketches...");

    // Create data structure to store combined sketches
    cluster_screen = new ArrayList<HashSet<Integer>>();
    for (int i = 0; i < Settings.NUM_CLUSTERS; i++){
      HashSet<Integer> tmp = new HashSet<Integer>();
      cluster_screen.add(tmp);
    }

    for (int x = 0; x < screen.genomeNames.length; x++)
    {
      HashSet<Integer> curr = screen.sketch_hash.get(x);

      // TODO - makee this a parameter
      HashSet<Integer> curr_shrunk = shrinkSet(curr, 1000);

      cluster_screen.get(assignments[x]).addAll(curr);

    }
  }

  HashSet<Integer> shrinkSet(HashSet<Integer> s, int shrinkSize) {

    ArrayList<Integer> set_list = new ArrayList<Integer>(s);

    HashSet<Integer> shrunk = new HashSet<Integer>();

    int size = s.size();

    for (int i = 0; i < shrinkSize; i++){
      int item = new Random().nextInt(size);
      int selected = set_list.get(item);
      shrunk.add(selected);
    }

    return shrunk;
  }

}
