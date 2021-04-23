package MessagesInDNA;

/* Contains functions
		patternMatching
		approxPatternMatching
		approxPatternCount
		hammingDistance
*/

import java.util.ArrayList;

public class PatternMatching{

	/* 	Function to find all occurences of a pattern in a segment of DNA
   		Input: String segment
   			   String pattern
   		Output: ArrayList of positions of each occurence of pattern 
   	*/
	public static ArrayList<Integer> patternMatching(String segment, String pattern){

		ArrayList<Integer> result = new ArrayList<>();

		for(int i = 0; i <= segment.length() - pattern.length(); i++){
			if(segment.substring(i,i+pattern.length()).equals(pattern))
				result.add(i);
		}

		return result;
	}

	/* 	Function to find all occurences of a pattern in a segment of DNA 
		given a threshold for the number of mismatches.
   		Input: String segment
   			   String pattern
   			   Integer d (max number of mismatches)
   		Output: ArrayList of positions of each occurence of pattern 
   	*/
   	public static ArrayList<Integer> approxPatternMatching(String segment, String pattern, Integer d){

   		int length = pattern.length();

   		ArrayList<Integer> result = new ArrayList<>();

		for(int i = 0; i <= segment.length() - pattern.length(); i++){

			int mismatches = 0;
			for(int j = 0; j < length; j++){
				if(segment.charAt(i+j) != pattern.charAt(j))
					mismatches++;
			}
			
			if(mismatches <= d)
				result.add(i);
		}
		return result;
   	}

   	/* 	Function to find the number of occurences of pattern in the given DNA
   		segment with a number of permitted mismatches d.
   		Input: String segment
   			   String pattern
   			   Integer d (max number of mismatches)
   		Output: ArrayList of positions of each occurence of pattern 
   	*/
   	public static Integer approxPatternCount(String segment, String pattern, Integer d){

		int count = 0;

		for(int i = 0; i <= segment.length() - pattern.length(); i++){

			String patternPrime = segment.substring(i,i+pattern.length());
			if(hammingDistance(pattern, patternPrime) <= d)
				count++;
		}
		return count;
   	}

	/* 	Function to find the hamming distance between two segments of DNA.
		The Hamming distance increases by 1 for every mismatch that occurs 
		between the nucleotides at position i in the two segments.
		Input: String segment1
			   String segment2
		Output: Integer hamming distance
	*/
	public static Integer hammingDistance(String s1, String s2){

		Integer hamming = 0;

		for(int i = 0; i < s1.length(); i++){
			if(s1.charAt(i) != s2.charAt(i))
				hamming++;
		}
		return hamming;
	}

}