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

// Heirachical Clustering Java
import com.apporiented.algorithm.clustering.*;
import com.apporiented.algorithm.clustering.visualization.*;

public class ParallelScreenerClusterNovel{

  // Concurrent List of reads to be processed
  ConcurrentLinkedQueue<String> reads_to_process;

  // Genome sketches
  ArrayList<HashSet<Integer>> sketch_hash;

  // Storing Cluster Information
  Cluster cluster;
  HashMap<String, ArrayList<String>> cluster_sketch_map;
  HashMap<String, Integer> genome_sketch_map;
  HashMap<String, Integer> cluster_height_map;
  HashMap<String, HashSet<Integer>> cluster_map;

  // Input set of reads to be processed
  ArrayList<String> reads;

  // Specified parameters
  int num_threads;
  int window;

  // Min Number of matches for read to be classified
  int threshold;

  // Source of these reads
  int readSet;

  // Total number of reads
  int totalReads;

  // Counts
  AtomicInteger insuf = new AtomicInteger(0);
  AtomicInteger ties = new AtomicInteger(0);

  // Keep track of reads
  AtomicInteger read_number;


  ParallelScreenerClusterNovel(ArrayList<HashSet<Integer>> sketch_hash, ArrayList<String> reads, int window, int readSet, int read_start, ClusterGenerator cg, HashMap<String, HashSet<Integer>> cluster_map){

    // Store parameters
    this.sketch_hash = sketch_hash;
    this.reads = reads;

    this.num_threads = Settings.NUM_THREADS;
    this.window = window;

    this.readSet = readSet;
    this.threshold = Settings.THRESHOLD;

    this.totalReads = 0;
    this.read_number = new AtomicInteger(read_start);


    // Cluster Information
    this.cluster = cg.cluster;
    this.cluster_sketch_map = cg.cluster_sketch_map;
    this.genome_sketch_map = cg.genome_sketch_map;
    this.cluster_height_map = cg.cluster_height_map;
    this.cluster_map = cluster_map;

    reads_to_process = new ConcurrentLinkedQueue<String>();
    // Populate queue
    int num_reads = reads.size();
    for (int r = 0; r < num_reads; r++){
      reads_to_process.add(reads.get(r));
      totalReads++;
    }
  }


  void run() throws Exception {

    MyThread[] threads = new MyThread[num_threads];
  		for(int i = 0; i < num_threads; i++)
  		{
  			threads[i] = new MyThread();
  			if(i == num_threads - 1)
  			{
  				threads[i].run();
  			} else {
  				threads[i].start();
  			}
  		}
  		for(int i = 0; i<num_threads-1; i++)
  		{
  			threads[i].join();
  		}
  }

  public class MyThread extends Thread {
    ;
    public void run() {

      while(!reads_to_process.isEmpty()) {

        String read = reads_to_process.poll();
        int curr_read = read_number.get();

        // Update read number
        read_number.incrementAndGet();

        ReadClassifierClusterNovel rc = new ReadClassifierClusterNovel(read, window, readSet, curr_read, cluster, cluster_sketch_map, genome_sketch_map, cluster_height_map, sketch_hash.size());
        rc.classifyRead(sketch_hash, cluster_map);

        // Update counts - TODO: this remains same for cluster
        insuf.addAndGet(rc.insufficient);
        ties.addAndGet(rc.tied);
      }

    }

  }
}
