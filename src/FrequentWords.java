package MessagesInDNA;

/* Contains functions:
		frequentWords
		frequentWordsWithMismatches
		frequentWordsWithMismatchesAndRC
		findNeighbors
*/

import java.util.HashMap;
import java.util.HashSet;
import java.util.stream.Collectors;
import java.util.Collections;

public class FrequentWords{


	/* 	Function to find the all k-mers of specified size k of a given DNA segment 
 	   	Function contains commented out alternate functionality to find the most 
 	   	frequent kmers in the given DNA segment
  	   	Input: String segment 
   			       Integer k 
  	   	Output: HashMap all k-mers, frequency pairs (or most frequent kmers)
	*/
	public static HashMap<String, Integer> frequentWords(String segment, Integer k){

		HashMap<String,Integer> frequencyTable = new HashMap<>();
		int maxFreq = 0;

		for(int i = 0; i < segment.length()-k; i++){
			/* This can be optimized by trimming/adding one character
			   at a time to the existing string instead of creating
			   a new one at each iteration */
			String pattern = segment.substring(i,i+k);

			frequencyTable.put(pattern,frequencyTable.getOrDefault(pattern,0)+1);

			if(frequencyTable.get(pattern) > maxFreq)
				maxFreq = frequencyTable.get(pattern);
		}

		return frequencyTable;
	}


	/* 	Function to find the most frequent k-mers of specified size k of a
		  given DNA segment with mismatch limit d. 
  	   	Input: String segment 
   			   Integer k 
   			   Integer d
  	   	Output: HashMap of most frequent kmers, count pairs
	*/
  	public static HashMap<String,Integer> frequentWordsWithMismatches(String segment, Integer k, Integer d){

  		HashMap<String,Integer> freqTable = new HashMap<>();
  		int maxFreq = 0;

  		for(int i = 0; i < segment.length()-k+1; i++){

  			String pattern = segment.substring(i,i+k);
  			HashSet<String> neighbors = findNeighbors(pattern,d);

  			for(String n : neighbors){
  				freqTable.put(n,freqTable.getOrDefault(n,0)+1);
  				if(freqTable.get(n) > maxFreq)
  					maxFreq = freqTable.get(n);
  			}
  		}

  		HashMap<String,Integer> mostFreqWords = new HashMap<>();
  		for(HashMap.Entry<String,Integer> item : freqTable.entrySet()){
  			if(item.getValue() == maxFreq)
  				mostFreqWords.put(item.getKey(),item.getValue());
  		}

  		return mostFreqWords;
  	}


  	/* 	Function to find the most frequent k-mers and their reverse
  		  complements of a given DNA segment with mismatch limit d. 
  	   	Input: String segment 
   			      Integer k 
   			      Integer d
  	   	Output: HashMap of most frequent kmers, count pairs
	*/
  	public static HashMap<String,Integer> frequentWordsWithMismatchesAndRC(String segment, Integer k, Integer d){

  		HashMap<String,Integer> freqTable = new HashMap<>();
  		int maxFreq = 0;

  		for(int i = 0; i < segment.length()-k+1; i++){

  			String pattern = segment.substring(i,i+k);
  			HashSet<String> neighbors = findNeighbors(pattern,d);

  			for(String n : neighbors){

  				String rcNeighbor = ReverseComplement.reverseComplement(n);

  				freqTable.put(n,freqTable.getOrDefault(n,0)+1);
  				freqTable.put(rcNeighbor,freqTable.getOrDefault(rcNeighbor,0)+1);

  				if(freqTable.get(n) > maxFreq)
  					maxFreq = freqTable.get(n);
  			}

  		}

  		HashMap<String,Integer> mostFreqWords = new HashMap<>();
  		for(HashMap.Entry<String,Integer> item : freqTable.entrySet()){
  			if(item.getValue() == maxFreq)
  				mostFreqWords.put(item.getKey(),item.getValue());
  		}

  		return mostFreqWords;
  	}


  	/* 	Function to find all variations of the given pattern with as many
		    as d mismatches.
		    Input: 	String pattern
				        Integer d (number of mismatches)
		    Output: Set of all neighbor patterns
	  */
  	public static HashSet<String> findNeighbors(String pattern, Integer d){

  		HashSet<String> nucleotides = new HashSet<>();
  		nucleotides.add("A");
  		nucleotides.add("T");
  		nucleotides.add("G");
		  nucleotides.add("C");

  		if(d == 0)
  			return (HashSet<String>) Collections.singleton(pattern);
  		if(pattern.length() == 1){
  			return nucleotides;
  		}

  		HashSet<String> neighborhood = new HashSet<>();

  		String firstSymbol = pattern.substring(0,1);
  		String suffix = pattern.substring(1);
  		HashSet<String> suffixNeighbors = findNeighbors(suffix,d);

  		for(String s : suffixNeighbors){
  			if(PatternMatching.hammingDistance(suffix,s) < d){
  				for(String n : nucleotides){
  					neighborhood.add(n + s);
  				}
  			}
  			else{
  				neighborhood.add(firstSymbol + s);
  			}
  		}
  		return neighborhood;
  	}


}