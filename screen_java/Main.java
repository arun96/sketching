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

// TODO - fix this eventually
public class Main {

  public static void main(String[] args) throws Exception {

    Settings.parseArgs(args);

    // Run the screen-generation and read classification processes
    run();

    // TODO - any post processing?

  }

  static void run() throws Exception {

    // TODO
    if (Settings.BAD_INPUT) {
      System.out.println("Invalid input parameters - please read the README!");
    } else {
      if (Settings.MINHASH){
        if (Settings.FIXED) {
          MinHashScreen screen = new MinHashScreen("fixed");
          ReadScreener rs = new ReadScreener(screen);
        } else {
          MinHashScreen screen = new MinHashScreen();
          ReadScreener rs = new ReadScreener(screen);
        }
      } else if (Settings.MINIMIZER){
        if (Settings.FIXED){
          MinimizerScreen screen = new MinimizerScreen("fixed");
          ReadScreener rs = new ReadScreener(screen);
        } else {
          MinimizerScreen screen = new MinimizerScreen();
          ReadScreener rs = new ReadScreener(screen);
        }

      } else if (Settings.UNIFORM) {
        UniformScreen screen = new UniformScreen();
        ReadScreener rs = new ReadScreener(screen);
      } else {
        System.out.println("Invalid input parameters - please read the README!");
      }
    }
  }
}
