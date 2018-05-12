package structures;

import shared.Tools;

public class SuperLongList {
	
	public SuperLongList(){
		this(100000);
	}
	
	public SuperLongList(int limit_){
		limit=limit_;
		array=new long[limit];
		list=new LongList();
	}
	
	public void add(long x){
		if(x<limit){array[(int)x]++;}
		else{list.add(x);}
		sum+=x;
		count++;
	}
	
	public void add(SuperLongList sllT){
		for(int i=0; i<sllT.array.length; i++){
			array[i]+=sllT.array[i];
		}
		list.append(sllT.list);
		count+=sllT.count;
		sum+=sllT.sum;
	}
	
	public double stdev(){
		final long div=Tools.max(1, count);
		double avg=sum/(double)div;
		double sumdev2=0;
		for(int i=0; i<array.length; i++){
			double dev=avg-i;
			double dev2=dev*dev;
			sumdev2+=(array[i]*dev2);
		}
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			double dev=avg-x;
			double dev2=dev*dev;
			sumdev2+=dev2;
		}
		return Math.sqrt(sumdev2/div);
	}
	
	/** Returns value such that percentile of values are below that value */
	public double percentileValueByCount(double percentile){
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentCount>=thresh){return i;}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev);
			currentSum+=x;
			currentCount++;
			if(currentCount>=thresh){return x;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}
	
	/** Returns value such that percentile of sum of values are below that value */
	public double percentileValueBySum(double percentile){
		long thresh=(long)(sum*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentSum>=thresh){return i;}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev);
			currentSum+=x;
			currentCount++;
			if(currentSum>=thresh){return x;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}
	
	public double percentileSumByCount(double percentile){
		long thresh=(long)(count*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentCount>=thresh){
				currentSum-=(x*i);
				currentCount-=i;
				while(currentCount<thresh){
					currentSum+=i;
					currentCount++;
				}
				return currentSum;
			}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev);
			currentSum+=x;
			currentCount++;
			if(currentCount>=thresh){return currentSum;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}
	
	public double percentileCountBySum(double percentile){
		long thresh=(long)(sum*percentile);
		long currentSum=0;
		long currentCount=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			currentSum+=(x*i);
			currentCount+=i;
			if(currentSum>=thresh){
				currentSum-=(x*i);
				currentCount-=i;
				while(currentSum<thresh){
					currentSum+=i;
					currentCount++;
				}
				return currentCount;
			}
		}
		long prev=-1;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			assert(x>=prev);
			currentSum+=x;
			currentCount++;
			if(currentSum>=thresh){return currentCount;}
			prev=x;
		}
		assert(false) : percentile+", "+count+", "+sum;
		return 0;
	}
	
	public long mode(){
		long maxCount=0;
		long maxValue=0;
		for(int i=0; i<array.length; i++){
			long x=array[i];
			if(x>maxCount){
				maxCount=x;
				maxValue=i;
			}
		}
		
		long prev=-1;
		long currentCount=0;
		for(int i=0; i<list.size; i++){
			long x=list.get(i);
			if(x==prev){
				currentCount++;
				if(currentCount>maxCount){
					maxCount=currentCount;
					maxValue=x;
				}
			}else{
				assert(x>prev) : "Needs to be sorted ascending.";
				prev=x;
				currentCount=1;
			}
		}
		return maxValue;
	}
	
	public void sort() {
		list.sort();
	}
	
	private long count;
	private long sum;
	
	final long[] array;
	final LongList list;
	final int limit;
	
}
