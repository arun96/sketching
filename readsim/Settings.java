/*
 * Settings for centromere assembler
 */
public class Settings {
	static int K = 13;
	static int MAX_KMER_MATCH_OFFSET = 50;
	static int MATCHES_NEEDED = 10;
	static int MAX_ITERS = 10;
	static int CENTROMER_COUNT = 70;
	static int MIN_UMER_COUNT = 20;
	static int MAX_UMER_COUNT = 165;
	static int MIN_READ_LENGTH = 500;
	static int MAX_ISLAND_SIZE = 987654321;
	static int MAX_SINGLE_READ = 10;
	static double UMER_PROPORTION = 0;
	static double CMER_PROPORTION = 0;
	static int MIN_OVERHANG = 1;
	static String UMER_SOURCE = "LR";
	static int MAX_NEIGHBORS = 10;
	
	static int NUM_THREADS = 1;

	static String READ_FILE = "";
	static String ASSEMBLY_FILE = "assembly.fa";
	static String LENGTHS_FILE = "lengths.txt";
    static String LAYOUT_FILE = "layout.txt";
    static String UMER_FILE = "umerlist.txt";

    static String PREFIX = "asm";
    static String EXT_READMER = ".readmer";
    static String EXT_FA = ".fa";
    static String EXT_LEN = ".len";
    static String EXT_LAY = ".lay";
    static String EXT_OVL = ".ovl";
    static String EXT_JOIN = ".join";
    static String EXT_UMER = ".umerlist";

    static Boolean VERBOSE = false;
	
	static void usage()
	{
		System.out.println();
		System.out.println("Centrotools");
		System.out.println("Usage: java -cp src Overlapper [args]");
		System.out.println("  Example: java -cp src Overlapper read_file=reads.fa");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  read_file (String) - a fasta file containing the reads");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  k                 [11] - kmer length used in detecting umers and cmers");
		System.out.println("  num_threads        [1] - number of threads to use for edge finding");
		System.out.println("  matches_needed     [6] - number of shared umers needed to join two reads or contigs");
		System.out.println("  max_offset       [100] - largest allowed change in relative kmer position when chaining umers");
		System.out.println("  max_island_size  [inf] - largest number of umers per read to combine into the same component");
		System.out.println("  max_single_read   [10] - largest number of times a umer can occur in a single read");
		System.out.println("  max_neighbors     [10] - number of neighbors to consider for each contig in a single iteration");
		System.out.println("  min_overhang     [500] - minimum amount of overhang to extend contig");
		System.out.println("  max_iters         [20] - maximum number of contig joining iterations");
		System.out.println("  min_umer_cnt      [30] - minimum coverage for a umer");
		System.out.println("  max_umer_cnt     [100] - maximum coverage for a umer");
		System.out.println("  cmer_cnt         [200] - minimum coverage for a cmer");
		System.out.println("  cmer_frac        [0.5] - fraction of cmers in read to be centromeric");
		System.out.println("  umer_frac       [0.03] - fraction of umers in read to be centromeric");
		System.out.println("  min_read_len     [500] - minimum read length to consider");
		System.out.println("  umer_source       [LR] - whether to get umers from long reads (LR), genome (G), or true mutations (T)");
		System.out.println("  assembly [assembly.fa] - filename for output assembly");
		System.out.println("  lengths  [lengths.txt] - filename for output contig lengths");
		System.out.println("  layout    [layout.txt] - filename for output contig layout");
        System.out.println("  prefix           [asm] - prefix for itermediate outputs");
        System.out.println("  -v                     - Activate verbose mode");
		System.out.println();
	}
	
	static boolean parseArgs(String[] args)
	{
		for(String s : args)
		{
            if(s.equalsIgnoreCase("-h") || s.equalsIgnoreCase("--help"))
            {
                usage();
                return false;
            }

            if(s.equalsIgnoreCase("-v") || s.equalsIgnoreCase("--verbose")) { VERBOSE = true; }

			int equalsIdx = s.indexOf('=');

			if(equalsIdx == -1)
			{
				continue;
			}
			else
			{
				String key = s.substring(0, equalsIdx);
				String val = s.substring(1 + equalsIdx);
                System.err.println("setting " + key + ": " + val);
                
				if(key.equalsIgnoreCase("read_file"))            { READ_FILE = val; }
				else if(key.equalsIgnoreCase("k"))               { K = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("num_threads"))               { NUM_THREADS = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("matches_needed"))  { MATCHES_NEEDED = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_offset"))      { MAX_KMER_MATCH_OFFSET = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_iters"))       { MAX_ITERS = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("min_umer_cnt"))    { MIN_UMER_COUNT = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_umer_cnt"))    { MAX_UMER_COUNT = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("cmer_cnt"))        { CENTROMER_COUNT = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("umer_frac"))       { UMER_PROPORTION = Double.parseDouble(val); }
				else if(key.equalsIgnoreCase("cmer_frac"))       { CMER_PROPORTION = Double.parseDouble(val); }
				else if(key.equalsIgnoreCase("min_read_len"))    { MIN_READ_LENGTH = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("min_overhang"))    { MIN_OVERHANG = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_island_size")) { MAX_ISLAND_SIZE = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_neighbors"))   { MAX_NEIGHBORS = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("max_single_read")) { MAX_SINGLE_READ = Integer.parseInt(val); }
				else if(key.equalsIgnoreCase("assembly"))        { ASSEMBLY_FILE = val; }
				else if(key.equalsIgnoreCase("lengths"))         { LENGTHS_FILE = val; }
				else if(key.equalsIgnoreCase("umer_source"))     { if(val.equals("G") || val.equals("T")) UMER_SOURCE = val; }

			}
		}
		
		if(READ_FILE.length() == 0)
		{
			usage();
			return false;
		}
		
		return true;
	}
}
