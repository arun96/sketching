import java.util.HashMap;

public class KmerCounter {
	HashMap<Long, Integer> counts;
	int k;
	long mask;

	KmerCounter(int k)
	{
		this.k = k;
		mask = (1L<<(2*k-2)) - 1;
		counts = new HashMap<Long, Integer>();
	}

	void addRead(String read)
	{
		long kmerVal = 0;
		for(int i = 0; i<read.length(); i++)
		{
			int cur = charToInt(read.charAt(i));
			kmerVal &= mask;
			kmerVal <<= 2;
			kmerVal |= cur;
			if(i >= k)
			{
				counts.put(kmerVal, 1 + counts.getOrDefault(kmerVal, 0));
			}
		}
	}

	int getFreq(long kmer)
	{
		return counts.getOrDefault(kmer, 0);
	}

	static int charToInt(char c)
	{
		if(c >= 'a' && c <= 'z') c += 'A' - 'a';
		if(c == 'A') return 0;
		else if(c == 'C') return 1;
		else if(c == 'G') return 2;
		return 3;
	}
	
	static char[] bases = new char[] {'A', 'C', 'G', 'T'};
	
	String kmerToString(long val)
	{
		char[] res = new char[k];
		for(int i = 0; i<k; i++)
		{
			res[k-1-i] = bases[(int)(val%4)];
			val /= 4;
		}
		return new String(res);
	}
	
	long stringToKmer(String s)
	{
		long res = 0;
		for(int i = 0; i<s.length(); i++)
		{
			res *= 4;
			res += charToInt(s.charAt(i));
		}
		return res;
	}

    void printStatus()
    {
      int umercnt = 0;
      int cmercnt = 0;
      int totalcnt = 0;

      for (long x : counts.keySet())
      {
        totalcnt++;

        int cnt = counts.get(x);

        if (cnt >= Settings.CENTROMER_COUNT) { cmercnt++; }
        else if ((cnt >= Settings.MIN_UMER_COUNT) && (cnt <= Settings.MAX_UMER_COUNT)) { umercnt++; }
      }

      System.err.println("Kmer counter mersloaded: " + totalcnt + " cmercnt: " + cmercnt + " umercnt: " + umercnt);
    }
}
