/*
 * Read simulator with reads coming from an exponential distribution
 */
import java.util.ArrayList;
import java.util.Random;
import java.io.File;
import java.io.PrintWriter;

public class ReadSimulator {
	String ref;
	double snpRate;
	double insRate;
	double delRate;
	int meanLength;
	int coverage;
	Random r;
	static int readId;
    boolean alwaysForward;

	ArrayList<String> names, reads, cigars;

	public static void main(String[] args) throws Exception
	{
		// Vanilla
		// String fn = args[0], ofn = args[1];
		// GenomeSimulator gs = new GenomeSimulator(fn);
		// ReadSimulator rs = new ReadSimulator(gs, 0.034, 0.033, 0.033, 10000, 10, false);
		// rs.printToFile(ofn);

		// adding command line arguments for error rates, genome length, and coverage
		String fn = args[0], ofn = args[1];

		String opt = "";

		// New Args
		double SNP = Double.parseDouble(args[2]);
		double insertion = Double.parseDouble(args[3]);
		double deletion = Double.parseDouble(args[4]);
		int meanLength = Integer.parseInt(args[5]);
		int coverage = Integer.parseInt(args[6]);
		if (args.length == 8) {
			opt = args[7];
		}

		GenomeSimulator gs = new GenomeSimulator(fn);
		ReadSimulator rs = new ReadSimulator(gs, SNP, insertion, deletion, meanLength, coverage, false, opt);
		rs.printToFile(ofn);
	}

	ReadSimulator(GenomeSimulator ref, double snpRate, double insRate, double delRate, int meanLength, int coverage, boolean alwaysForward, String opt)
	{
		r = new Random(654321);
		this.ref = new String(ref.s);
		this.snpRate = snpRate;
		this.delRate = delRate;
		this.insRate = insRate;
		this.meanLength = meanLength;
		this.coverage = coverage;
        this.alwaysForward = alwaysForward;

		names = new ArrayList<String>();
		reads = new ArrayList<String>();
		cigars = new ArrayList<String>();

		long totLength = 0;
		long lengthNeeded = (long) ref.s.length * coverage;
		while(totLength < lengthNeeded)
		{
			String[] curRead = makeRead(opt);
			names.add(curRead[0]);
			reads.add(curRead[1]);
			cigars.add(curRead[2]);
			totLength += curRead[1].length() * 1.0 / (1 - delRate + insRate);
		}
	}

    void printToFile(String filename) throws Exception
    {
      PrintWriter out = new PrintWriter(new File(filename));

      for (int i = 0; i < reads.size(); i++)
      {
        out.println(">" + names.get(i));
        out.println(reads.get(i));
      }

      out.close();
    }

	String[] makeRead(String opt)
	{
		int length = generateReadLength(opt);
		if(length > ref.length())
		{
			length = ref.length();
		}
		int start = r.nextInt(ref.length() - length + 1);
		int end = start + length;
		char strand = (r.nextInt(2) == 0) ? '+' : '-';

        if (alwaysForward) { strand = '+'; }

		String[] readSeq = simulateRead(start, end, strand);
		String name = "read" + readId + "_" + start + "_" + end + "_" + readSeq[0].length() + "_" + (strand == '+' ? 'F' : 'R');
		readId++;
		return new String[] {name, readSeq[0], readSeq[1]};
	}

	int generateReadLength(String opt)
	{

		if (opt.equals("E")){
			// Original Exponential
			return (int)(.5 + Math.log(1-r.nextDouble())*(-meanLength));
		} else if (opt.equals("EM")){
			// Exponential with Min Length
			return (int)(Math.max((.5 + Math.log(1-r.nextDouble())*(-meanLength)), (meanLength/2.0)));
		} else if (opt.equals("EL")){
			// Looping Exponential
			double readLen = 0.00;
			while (readLen < meanLength/2.0)
			{
				readLen = (.5 + Math.log(1-r.nextDouble())*(-meanLength));
			}
			return (int) (readLen);
		} else if (opt.equals("XL")){
			// Exact Read Lengths
			return meanLength;

		} else {
			//  DEFAULT - Normal Distribution
			return (int)(r.nextGaussian()*(Math.sqrt(meanLength))+meanLength);
		}


	}

	static String revComp(String s)
	{
		int n = s.length();
		char[] res = new char[n];
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(n - 1 - i);
			if(c == 'a' || c == 'A') c += 'T' - 'A';
			else if(c == 'c' || c == 'C') c += 'G' - 'C';
			else if(c == 'g' || c == 'G') c += 'C' - 'G';
			else c += 'A' - 'T';
			res[i] = c;
		}
		return new String(res);
	}

	String getRefSeq(int start, int end, char strand)
	{
		String res = ref.substring(start, end);
		if(strand == '+')
		{
			return res;
		}
		else
		{
			return revComp(res);
		}
	}

	char mutate(char c)
	{
		char[] options = (c >= 'A' && c <= 'Z') ? new char[] {'A', 'C', 'G', 'T'} : new char[] {'a', 'c', 'g', 't'};
		char res = options[r.nextInt(4)];
		while(res == c)
		{
			res = options[r.nextInt(4)];
		}
		return res;
	}

	char randomChar()
	{
		char[] options = new char[] {'A', 'C', 'G', 'T'};
		return options[r.nextInt(4)];
	}

	String[] simulateRead(int start, int end, char strand)
	{
		String s = getRefSeq(start, end, strand);
		int idx = 0;
		StringBuilder sb = new StringBuilder("");
		StringBuilder cigar = new StringBuilder("");
		while(idx < s.length())
		{
			char cur = s.charAt(idx);
			double type = r.nextDouble();
			if(type < snpRate)
			{
				cigar.append('X');
				sb.append(mutate(cur));
				idx++;
			}
			else if(type < snpRate + insRate)
			{
				cigar.append('I');
				sb.append(randomChar());
			}
			else if(type < snpRate + insRate + delRate)
			{
				cigar.append('D');
				idx++;
			}
			else
			{
				cigar.append('M');
				sb.append(cur);
				idx++;
			}
		}
		return new String[] {sb.toString(), cigar.toString()};
	}
}
