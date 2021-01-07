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

public class ParallelScreenerNovel{

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

  // Total number of reads
  int totalReads;

  // The readset being classified
  int readSet;

  // Counts
  AtomicInteger insuf = new AtomicInteger(0);
  AtomicInteger ties = new AtomicInteger(0);

  // Keep track of reads
  AtomicInteger read_number;


  ParallelScreenerNovel(ArrayList<HashSet<Integer>> sketch_hash, int readSet, ArrayList<String> reads, int window, int read_start){

    // Store parameters
    this.sketch_hash = sketch_hash;
    this.reads = reads;

    this.num_threads = Settings.NUM_THREADS;
    this.window = window;

    this.threshold = Settings.THRESHOLD;

    this.totalReads = 0;
    this.read_number = new AtomicInteger(read_start);
    this.readSet = readSet;

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

        ReadClassifierNovel rc = new ReadClassifierNovel(sketch_hash, readSet, read, window, curr_read);
        rc.classifyRead();

        // Update counts
        insuf.addAndGet(rc.insufficient);
        ties.addAndGet(rc.tied);
      }

    }

  }
}