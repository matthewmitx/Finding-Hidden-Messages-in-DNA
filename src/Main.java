import MessagesInDNA.FrequentWords;
import MessagesInDNA.MotifEnumeration;

import java.util.*;
import java.io.*;

class Main{

	/* Main used for frequent words with mismatches program */
	public static void main(String[] args){

		Scanner scanner = new Scanner(System.in);

		Integer k = scanner.nextInt();
		scanner.nextLine();

		List<String> dna = new ArrayList<>();

		while(scanner.hasNext()){
			dna.add(scanner.nextLine());
		}

		System.out.println(MotifEnumeration.medianString(k,dna));

	}

}