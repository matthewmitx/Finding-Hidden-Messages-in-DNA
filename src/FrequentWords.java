package MessagesInDNA;

/* Contains functions:
		frequentWords
		hammingDistance
*/

import java.util.HashMap;

public class FrequentWords{


	/* 	Function to find the all k-mers of specified size k of a given DNA segment 
 	   	Function contains commented out alternate functionality to find the most 
 	   	frequent kmers in the given DNA segment
  	   	Input: String segment 
   			   Integer k 
  	   	Output: HashMap all k-mers, frequency pairs (or most frequent kmers)
	*/
	public static HashMap<String, Integer> frequentWords(String segment, int k){

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

		/* Modification to find the most frequent kmers */
		/* HashMap<String,Integer> result = new HashMap<>();
		for(HashMap.Entry<String,Integer> entry : frequencyTable.entrySet()){
			
    		if(entry.getValue() == maxFreq){
    			result.add(entry.getKey(),entry.getValue());
    		}
		} */
	}

}