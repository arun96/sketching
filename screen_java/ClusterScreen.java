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

      // TODO: Downsample before we add
      
      cluster_screen.get(assignments[x]).addAll(curr);

    }
  }

}
