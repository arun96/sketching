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

// Main function
public class Main {

  public static void main(String[] args) throws Exception {

    Settings.parseArgs(args);

    // System.out.println(Settings.READ_LENGTH + " " + Settings.READ_ERROR + " " + Settings.TARGET_MATCHES);

    if (Settings.SCREEN_ONLY) {
      run_screen();
    }

    else {
      // Run the screen-generation and read classification processes
      run();
    }

  }

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
        UniformScreen screen = new UniformScreen();

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
        UniformScreen screen = new UniformScreen();
        SaveScreen SS = new SaveScreen(screen);
      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }
}
