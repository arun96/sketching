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
import java.util.Collection;

public class ExhaustiveScreen extends Screen {

  public static void main(String[] args) throws Exception
  {
  }

  // Screen Generator for MinHash-based screen
  ExhaustiveScreen() throws Exception
  {
    System.out.println("Generating Exhaustive Screen, using " + getHashName(Settings.HASH_TYPE) + " hash function.");

    // Store variables
    this.targetMatches = Settings.TARGET_MATCHES;
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.GENOME_FOLDER;
    this.k = Settings.K;
    this.genomeNames = Settings.GENOMES;
    this.numGenomes = genomeNames.length;
    this.window = 0;

    // Array for the genomes
    String[] genomes = new String[numGenomes];

    // Read in the genomes, save length and the genome
    for (int i = 0; i < numGenomes; i++)
    {
      String gnm = getGenome(genomeFolder + genomeNames[i]);
      genomes[i] = gnm;
    }

    // Store the sketch
    sketch_hash = new ArrayList<HashSet<Integer>>();

    // For each genome
    for (int x = 0; x < numGenomes; x++)
    {
      // Get all hashes in the genome
      sketch_hash.add(getAllHashes(genomes[x], k));
    }
  }

  // Option for pre-generated screens
  ExhaustiveScreen(String[] genomes, ArrayList<HashSet<Integer>> sketch) throws Exception {
    System.out.println("Creating Exhaustive Screen...");
    this.targetMatches = Settings.TARGET_MATCHES;
    this.readLen = Settings.READ_LENGTH;
    this.readErr = Settings.READ_ERROR;
    this.genomeFolder = Settings.SCREEN_LOCATION;
    this.k = Settings.K;
    this.genomeNames = genomes;
    this.numGenomes = genomeNames.length;
    this.window = 0;
    sketch_hash = sketch;
  }

}
