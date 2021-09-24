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

// Main function
public class Main {

  public static void main(String[] args) throws Exception {

    // Parse inputs and populate settings variables
    Settings.parseArgs(args);

    // If the genomes and reads are paired, we can  report accuracy/know the ground truth
    if (Settings.MATCHED_READS_GENOMES) {
      System.out.println("Running experiment with matched genomes and reads - accuracy will be computed and output.");
    } else {
      System.out.println("Running experiment unmatched genomes and reads - read classification results can be found in log files.");
    }

    // If Read Logging is enabled
    if (Settings.READ_LOGGING){
      System.out.println("Read Logging is enabled - logs will be saved here: " + Settings.READ_LOCATION);
    }

    // Runs the read screening
    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
      return;
    } else {
      if (Settings.SCREEN_ONLY) {
        System.out.println("Screen-Generation mode - screens will be saved here: " + Settings.SCREEN_LOCATION);
        run_screen();
      } else if (Settings.LOAD_SCREEN) {
        System.out.println("Loading pre-generated screen from: " + Settings.SCREEN_LOCATION + ". Please make sure input parameters match those of generated screens (found in params.txt).");
        run_load();
      } else if (Settings.CLUSTER_BASED) {
        run_cluster();
      } else {
        // Run the screen-generation and read classification processes
        run();
      }
    }
  }

  // ---- Vanilla screen-generation and read-screening ----
  static void run() throws Exception {

    Screen screen = null;

    // MinHash
    if (Settings.MINHASH){
      screen = new MinHashScreen();

    // Minimizer
    } else if (Settings.MINIMIZER){
      screen = new MinimizerScreen();

    // Uniform
    } else if (Settings.UNIFORM) {
      screen = new UniformScreen();

    // Exhaustive
    } else if (Settings.EXHAUSTIVE) {
      screen = new ExhaustiveScreen();

    } else {
      System.out.println("Invalid input parameters - please read the README!");
      return;
    }

    if (Settings.SAVE_SCREEN){
      SaveScreen SS = new SaveScreen(screen);
    }

    if (Settings.MATCHED_READS_GENOMES) {
      ReadScreener rs = new ReadScreener(screen);
    } else {
      ReadScreenerNovel rs = new ReadScreenerNovel(screen);
    }
  }

  // ---- Save screens, does not do read classification ----
  static void run_screen() throws Exception {

    System.out.println("Generating and saving screen!");

    Screen screen = null;

    if (Settings.MINHASH){
      screen = new MinHashScreen();

    } else if (Settings.MINIMIZER){
      screen = new MinimizerScreen();

    } else if (Settings.UNIFORM) {
      screen = new UniformScreen();

    } else if (Settings.EXHAUSTIVE) {
      screen = new ExhaustiveScreen();

    } else {
      System.out.println("Invalid input parameters - please read the README!");
      return;
    }

    SaveScreen SS = new SaveScreen(screen);
  }

  // ---- Load a pre-generated screen, and then run regular read screening. ----
  static void run_load() throws Exception {

    // Load saved screen
    LoadScreen ls = new LoadScreen();

    Screen screen = null;

    if (Settings.MINHASH) {
      screen = new MinHashScreen(ls.genomeNames, ls.sketch, ls.weights);

    } else if (Settings.MINIMIZER){
      screen = new MinimizerScreen(ls.genomeNames, ls.sketch);

    } else if (Settings.UNIFORM) {
      screen = new UniformScreen(ls.genomeNames, ls.sketch);

    } else if (Settings.EXHAUSTIVE) {
      screen = new ExhaustiveScreen(ls.genomeNames, ls.sketch);

    } else {
      System.out.println("Invalid input parameters - please read the README!");
      return;
    }

    if (Settings.MATCHED_READS_GENOMES) {
      ReadScreener rs = new ReadScreener(screen);
    } else {
      ReadScreenerNovel rs = new ReadScreenerNovel(screen);
    }
  }

  // ---- Run the clustering based approach ----
  static void run_cluster() throws Exception {

    // Cluster generation
    System.out.println("Clustering genomes...");
    ClusterGenerator cg = new ClusterGenerator();
    System.out.println("Clusters generated:");

    // DEBUGGING
    cg.cluster.toConsole(2);
    // System.out.println(cg.genome_sketch_map.toString());

    System.out.println("Starting Screen Generation and Classification...");

    Screen screen = null;

    if (Settings.MINHASH){
      screen = new MinHashScreen();

    } else if (Settings.MINIMIZER){
      screen = new MinimizerScreen();

    } else if (Settings.UNIFORM) {
      screen = new UniformScreen();

    } else if (Settings.EXHAUSTIVE) {
      screen = new ExhaustiveScreen();

    } else {
      System.out.println("Invalid input parameters - please read the README!");
      return;
    }

    if (Settings.SAVE_SCREEN){
      SaveScreen SS = new SaveScreen(screen);
    }

    if (Settings.MATCHED_READS_GENOMES) {
      ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
    } else {
      ReadScreenerClusterNovel rs = new ReadScreenerClusterNovel(screen, cg);
    }
  }
}
