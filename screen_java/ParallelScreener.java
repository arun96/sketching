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

public class ParallelScreener{

  // Concurrent List of reads to be processed
  ConcurrentLinkedQueue<String> reads_to_process;

  // Genome sketches
  ArrayList<HashSet<Integer>> sketch_hash;

  // Input set of reads to be processed
  ArrayList<String> reads;

  // Specified parameters
  int num_threads;
  int window;

  // Min Number of matches for read to be classified
  int threshold;

  // Source of these reads
  int source;

  // Total number of reads
  int totalReads;

  // Counts
  AtomicInteger correct = new AtomicInteger(0);
  AtomicInteger mis = new AtomicInteger(0);
  AtomicInteger insuf = new AtomicInteger(0);
  AtomicInteger ties = new AtomicInteger(0);

  // Keep track of reads
  AtomicInteger read_number;


  ParallelScreener(ArrayList<HashSet<Integer>> sketch_hash, ArrayList<String> reads, int window, int source, int read_start){

    // Store parameters
    this.sketch_hash = sketch_hash;
    this.reads = reads;

    this.num_threads = Settings.NUM_THREADS;
    this.window = window;

    this.source = source;
    this.threshold = Settings.THRESHOLD;

    this.totalReads = 0;
    this.read_number = new AtomicInteger(read_start);

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

        ReadClassifier rc = new ReadClassifier(sketch_hash, read, window, source, curr_read);
        rc.classifyRead();

        // Update counts
        correct.addAndGet(rc.correct);
        mis.addAndGet(rc.incorrect);
        insuf.addAndGet(rc.insufficient);
        ties.addAndGet(rc.tied);
      }

    }

  }






}
