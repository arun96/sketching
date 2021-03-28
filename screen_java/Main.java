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

    Settings.parseArgs(args);

    if (Settings.MATCHED_READS_GENOMES) {
      System.out.println("Running experiment with matched genomes and reads - accuracy will be computed and output.");
    } else {
      System.out.println("Running experiment unmatched genomes and reads - read classification results can be found in log files.");
    }

    if (Settings.READ_LOGGING){
      System.out.println("Read Logging is enabled - logs will be saved here: " + Settings.READ_LOCATION);
    }

    // System.out.println(Settings.READ_LENGTH + " " + Settings.READ_ERROR + " " + Settings.TARGET_MATCHES);

    if (Settings.SCREEN_ONLY) {
      System.out.println("Screen-Generation mode - screens will be saved here: " + Settings.SCREEN_LOCATION);
      run_screen();
    } else if (Settings.LOAD_SCREEN) {
      System.out.println("Loading pre-generated screen from: " + Settings.SCREEN_LOCATION);
      run_load();
    } else if (Settings.CLUSTER_BASED) {
      run_cluster();
    } else {
      // Run the screen-generation and read classification processes
      run();
    }

  }

  // Vanilla screen-generation and read-screening
  static void run() throws Exception {
    // Runs the read screening
    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
    } else {
      if (Settings.MINHASH){
        if (Settings.FIXED) {
          MinHashScreen screen = new MinHashScreen("fixed");

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        } else {
          MinHashScreen screen = new MinHashScreen();

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        }
      } else if (Settings.MINIMIZER){
        if (Settings.FIXED){
          MinimizerScreen screen = new MinimizerScreen("fixed");

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        } else {
          MinimizerScreen screen = new MinimizerScreen();

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        }

      } else if (Settings.UNIFORM) {
        if (Settings.FIXED){
          UniformScreen screen = new UniformScreen("fixed");

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }
        } else {
          UniformScreen screen = new UniformScreen();

          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreener rs = new ReadScreener(screen);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }
        }

      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }

  // Save screens (either as a whole, or individually)
  static void run_screen() throws Exception {

    System.out.println("Generating and saving screen!");

    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
    } else {
      if (Settings.MINHASH){
        if (Settings.FIXED) {
          MinHashScreen screen = new MinHashScreen("fixed");
          SaveScreen SS = new SaveScreen(screen);
        } else {
          MinHashScreen screen = new MinHashScreen();
          SaveScreen SS = new SaveScreen(screen);
        }
      } else if (Settings.MINIMIZER){
        if (Settings.FIXED){
          MinimizerScreen screen = new MinimizerScreen("fixed");
          SaveScreen SS = new SaveScreen(screen);
        } else {
          MinimizerScreen screen = new MinimizerScreen();
          SaveScreen SS = new SaveScreen(screen);
        }

      } else if (Settings.UNIFORM) {
        if (Settings.FIXED){
          UniformScreen screen = new UniformScreen("fixed");
          SaveScreen SS = new SaveScreen(screen);
        } else {
          UniformScreen screen = new UniformScreen();
          SaveScreen SS = new SaveScreen(screen);
        }
      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }

  // Load a pre-generated screen, and then run regular read screening.
  static void run_load() throws Exception {
    // Runs the read screening
    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
    } else {
      LoadScreen ls = new LoadScreen();
      if (Settings.MINHASH) {
        MinHashScreen screen = new MinHashScreen(ls.genomeNames, ls.sketch);

        if (Settings.MATCHED_READS_GENOMES) {
          ReadScreener rs = new ReadScreener(screen);
        } else {
          ReadScreenerNovel rs = new ReadScreenerNovel(screen);
        }

      } else if (Settings.MINIMIZER){
        MinimizerScreen screen = new MinimizerScreen(ls.genomeNames, ls.sketch);

        if (Settings.MATCHED_READS_GENOMES) {
          ReadScreener rs = new ReadScreener(screen);
        } else {
          ReadScreenerNovel rs = new ReadScreenerNovel(screen);
        }

      } else if (Settings.UNIFORM) {
        UniformScreen screen = new UniformScreen(ls.genomeNames, ls.sketch);

        if (Settings.MATCHED_READS_GENOMES) {
          ReadScreener rs = new ReadScreener(screen);
        } else {
          ReadScreenerNovel rs = new ReadScreenerNovel(screen);
        }


      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }

  static void run_cluster() throws Exception {
    System.out.println("Clustering genomes...");
    ClusterGenerator cg = new ClusterGenerator();
    System.out.println("Clusters generated:");

    // DEBUGGING
    cg.cluster.toConsole(2);
    //System.out.println(cg.cluster_sketch_map.toString());
    // System.out.println(cg.genome_sketch_map.toString());

    System.out.println("Starting Screen Generation and Classification...");

    // Runs the read screening

    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
    } else {
      if (Settings.MINHASH){
        if (Settings.FIXED) {
          MinHashScreen screen = new MinHashScreen("fixed");


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        } else {
          MinHashScreen screen = new MinHashScreen();


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        }
      } else if (Settings.MINIMIZER){
        if (Settings.FIXED){
          MinimizerScreen screen = new MinimizerScreen("fixed");


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        } else {
          MinimizerScreen screen = new MinimizerScreen();


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }

        }

      } else if (Settings.UNIFORM) {
        if (Settings.FIXED){
          UniformScreen screen = new UniformScreen("fixed");


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }
        } else {
          UniformScreen screen = new UniformScreen();


          if (Settings.MATCHED_READS_GENOMES) {
            ReadScreenerCluster rs = new ReadScreenerCluster(screen, cg);
          } else {
            ReadScreenerNovel rs = new ReadScreenerNovel(screen);
          }
        }
        // STOP

      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }
}
