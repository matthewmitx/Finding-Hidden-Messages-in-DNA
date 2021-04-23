import MessagesInDNA.FrequentWords;
import MessagesInDNA.MinimumSkew;
import MessagesInDNA.FindClumps;

import java.util.*;
import java.io.*;

class Main{

	/* Main used for frequent words with mismatches program */
	public static void main(String[] args){

		HashSet<String> n = FrequentWords.findNeighbors("CCAGTCAATG", 1);
		System.out.println(n.size());
		System.out.println(n.toString());

		/*Scanner scanner = new Scanner(System.in);
		String genome = scanner.next();

		ArrayList<Integer> mins = MinimumSkew.minimumSkew(genome);

		Integer minimum = mins.get(0);

		System.out.println(minimum);

		String oriCentered = genome.substring(minimum-250,minimum+250);
		String oriLeft = genome.substring(minimum-500,minimum);
		String oriRight = genome.substring(minimum,minimum+500);

		HashMap<String,Integer> dnaABoxesCentered = FrequentWords.frequentWordsWithMismatchesAndRC(oriCentered,9,1);
		System.out.println("Ori centered: " + dnaABoxesCentered.size());
		for(HashMap.Entry<String,Integer> item : dnaABoxesCentered.entrySet()){
			System.out.println(item.getKey() + ": " + item.getValue());
		}

		HashMap<String,Integer> dnaABoxesLeft = FrequentWords.frequentWordsWithMismatchesAndRC(oriLeft,9,1);
		System.out.println("Ori left: " + dnaABoxesLeft.size());
		for(HashMap.Entry<String,Integer> item : dnaABoxesLeft.entrySet()){
			System.out.println(item.getKey() + ": " + item.getValue());
		}

		HashMap<String,Integer> dnaABoxesRight = FrequentWords.frequentWordsWithMismatchesAndRC(oriRight,9,1);
		System.out.println("Ori right: " + dnaABoxesRight.size());
		for(HashMap.Entry<String,Integer> item : dnaABoxesRight.entrySet()){
			System.out.println(item.getKey() + ": " + item.getValue());
		}*/

	}

}