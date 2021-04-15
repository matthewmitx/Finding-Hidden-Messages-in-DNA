/* Program to find all occurences of a pattern in a segment of DNA
   Input: String segment
   		  String pattern
   Output: ArrayList of positions of each occurence of pattern */

import java.util.ArrayList;

class PatternMatching{

	public static ArrayList<Integer> patternMatching(String segment, String pattern){

		ArrayList<Integer> result = new ArrayList<>();

		for(int i = 0; i < segment.length() - pattern.length(); i++){
			if(segment.substring(i,i+pattern.length()).equals(pattern))
				result.add(i);
		}

		return result;
	}

}