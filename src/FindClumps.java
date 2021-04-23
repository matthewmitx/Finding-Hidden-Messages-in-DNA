package MessagesInDNA;

/* Contains functions
		findClumps
		findClumpsWithMismatchesAndRC
*/

import java.util.HashMap;
import java.util.HashSet;

public class FindClumps{

	/* 	Function to find all l sized clumps in a DNA segment that have any k-mer that 
   		repeats at least t times.
   		Input: String segment - full DNA segment to traverse
   		  	int k - size of kmer pattern
   		  	int l - size of clump in which to find repeating patterns
   		  	int t - minimum number of repititions for qualifying kmer pattern
   		Output: All qualifying kmers in segment
   		   		(size k patterns in repeating at least t times in a window of size l)
	*/
	public static HashSet<String> findClumps(String segment, int k, int l, int t){

		HashSet<String> patterns = new HashSet<>();

		int n = segment.length();

		/* Implementation that updates the set of qualifying kmer patters by 
		   removing on occurence of a single pattern and adding one occurence
		   of a single patter on each iteration without repeated traversal
		   Time Complexity: O()
		*/
		String window = segment.substring(0,l);
		HashMap<String,Integer> freqMap = FrequentWords.frequentWords(window,k);

		for(HashMap.Entry<String,Integer> entry : freqMap.entrySet()){
			if(entry.getValue() > t){
				patterns.add(entry.getKey());
			}
		}

		for(int i = 1; i < n-l; i++){
			String oldKmer = segment.substring(i-1,i+k-1);
			String newKmer = segment.substring(i+l-k,i+l);

			if(freqMap.containsKey(oldKmer)){
				freqMap.replace(oldKmer,freqMap.get(oldKmer)-1);
			}
			freqMap.put(newKmer,freqMap.getOrDefault(newKmer,0)+1);
			if(freqMap.get(newKmer) >= t && !patterns.contains(newKmer)){
				patterns.add(newKmer);
			}
		}
		return patterns;
	}


	/* !!!!!!!!!! THIS HAS NOT YET BEEN TESTED !!!!!!!!!! */

	/* 	Function to find all l sized clumps in a DNA segment that have any k-mer that 
   		repeats at least t times. K-mers are allowed to have d mismatches.
   		Input: String segment - full DNA segment to traverse
   		  	int k - size of kmer pattern
   		  	int d - mismatch limit for each kmer
   		  	int l - size of clump in which to find repeating patterns
   		  	int t - minimum number of repititions for qualifying kmer pattern
   		Output: All qualifying kmers in segment
   		   		(size k patterns in repeating at least t times in a window of size l)
	*/
	public static HashSet<String> findClumps(String segment, int k, int d, int l, int t){

		HashSet<String> patterns = new HashSet<>();

		int n = segment.length();

		String window = segment.substring(0,l);
		HashMap<String,Integer> freqMap = FrequentWords.frequentWordsWithMismatchesAndRC(window,k,d);

		for(HashMap.Entry<String,Integer> entry : freqMap.entrySet()){
			if(entry.getValue() > t){
				patterns.add(entry.getKey());
			}
		}

		for(int i = 1; i < n-l; i++){
			String oldKmer = segment.substring(i-1,i+k-1);
			String newKmer = segment.substring(i+l-k,i+l);

			HashSet<String> oldNeighbors = FrequentWords.findNeighbors(oldKmer,d);
			HashSet<String> newNeighbors = FrequentWords.findNeighbors(newKmer,d);

			for(String oldN : oldNeighbors){
				if(freqMap.containsKey(oldN)){
					freqMap.replace(oldN,freqMap.get(oldN)-1);
				}
			}

			for(String newN : newNeighbors){
				freqMap.put(newN,freqMap.getOrDefault(newN,0)+1);
				if(freqMap.get(newN) >= t && !patterns.contains(newN)){
					patterns.add(newN);
				}
			}
		}
		return patterns;
	}
	

}