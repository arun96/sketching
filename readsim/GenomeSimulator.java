/*
 * Simple centromeric genome simulator
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.Scanner;
import java.util.TreeSet;

public class GenomeSimulator {
	
	int flankingSequence;
	int repeatCount;
	int repeatLength;
	double mutationRate;
	TreeSet<Integer> mutations;
	char[] s;
	Random r;
	
	public static void main(String[] args) throws Exception
	{
		String fn = args[0];
		String ofn = args[1];
		GenomeSimulator gs = new GenomeSimulator(fn);
		gs.printToFile(ofn);
	}
	
	GenomeSimulator reverseComplement()
	{
		GenomeSimulator res = new GenomeSimulator();
		int n = s.length;
		res.s = ReadSimulator.revComp(new String(s)).toCharArray();
		res.mutations = new TreeSet<Integer>();
		for(int m : mutations)
		{
			res.mutations.add(n - 1 - m);
		}
		res.flankingSequence = this.flankingSequence;
		res.repeatCount = this.repeatCount;
		res.repeatLength = this.repeatLength;
		res.mutationRate = this.mutationRate;
		return res;
	}
	
	GenomeSimulator()
	{
		
	}
	
	GenomeSimulator(double mutationRate, int repeatLength, int repeatCount, int flankingSequence)
	{
		this.mutationRate = mutationRate;
		this.repeatLength = repeatLength;
		this.repeatCount = repeatCount;
		this.flankingSequence = flankingSequence;
		
		r = new Random(123);
		mutations = new TreeSet<Integer>();
		s = new char[flankingSequence*2 + repeatLength * repeatCount];
		int idx = 0;
		for(int i = 0; i<flankingSequence; i++)
		{
			s[idx++] = randomBase();
		}
		char[] repeatSeq = new char[repeatLength];
		for(int i = 0; i<repeatLength; i++)
		{
			repeatSeq[i] = randomBase();
		}
		for(int i = 0; i<repeatCount; i++)
		{
			for(int j = 0; j<repeatLength; j++)
			{
				s[idx] = repeatSeq[j];
				if(r.nextDouble() < mutationRate)
				{
					s[idx] = flipBase(s[idx]);
					mutations.add(idx);
				}
				idx++;
			}
		}
		for(int i = 0; i<flankingSequence; i++)
		{
			s[idx++] = randomBase();
		}
	}
	
	GenomeSimulator(String fn) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		StringBuilder sb = new StringBuilder("");
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith(">"))
			{
				continue;
			}
			sb.append(line);
		}
		String s = sb.toString();
		this.s = s.toUpperCase().toCharArray();
		this.mutations = new TreeSet<Integer>();
		for(int i = 0; i<s.length(); i++)
		{
			if(s.charAt(i) >= 'a' && s.charAt(i) <= 'z')
			{
				this.mutations.add(i);
			}
		}
		input.close();
	}
	
	void printToFile(String fn) throws Exception
	{
		PrintWriter out = new PrintWriter(new File(fn));
		out.println(">genome");
		for(int i = 0; i<s.length; i++)
		{
			char toPrint = s[i];
			if(mutations.contains(i))
			{
				toPrint = (char)(toPrint + 'a' - 'A');
			}
			out.print(toPrint);
			if(i%80 == 79 || i == s.length - 1)
			{
				out.println();
			}
		}
		out.close();
	}
	
	/*
	 * Prints the unique kmers in the genome to a file
	 */
	void printUmersToFile(String fn) throws Exception
	{
		KmerCounter kmc = new KmerCounter(Settings.K);
		String genomeString = new String(s);
		kmc.addRead(new String(s));
		
		ArrayList<String> umers = new ArrayList<String>();
		for(int i = 0; i + Settings.K <= repeatCount * repeatLength; i++)
		{
			String cur = genomeString.substring(i + flankingSequence, i + flankingSequence + Settings.K);
			int count = kmc.getFreq(kmc.stringToKmer(cur));
			if(count == 1)
			{
				umers.add(cur);
			}
		}
		Collections.sort(umers);
		
		PrintWriter out = new PrintWriter(new File(fn));
		for(String s : umers)
		{
			out.println(s+"\t"+1);
		}
		out.close();
	}
	
	static char[] bases = new char[] {'A', 'C', 'G', 'T'};
	char randomBase()
	{
		int val = r.nextInt(4);
		return bases[val];
	}
	char flipBase(char c)
	{
		if(c == 'A') return 'T';
		else if(c == 'C') return 'G';
		else if(c == 'G') return 'C';
		else return 'A';
	}
}
