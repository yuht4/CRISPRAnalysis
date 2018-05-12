package structures;

public class LongPair implements Comparable<LongPair>{

	@Override
	public int compareTo(LongPair other) {
		if(a!=other.a){return a>other.a ? 1 : -1;}
		return b>other.b ? 1 : b<other.b ? -1 : 0;
	}
	
	public long a, b;
	
}
