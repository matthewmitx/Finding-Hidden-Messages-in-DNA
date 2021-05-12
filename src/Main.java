import MessagesInDNA.FrequentWords;
import MessagesInDNA.MotifEnumeration;
import MessagesInDNA.Helper;

import java.util.*;
import java.io.*;

class Main{

	/* Main used for frequent words with mismatches program */
	public static void main(String[] args){

		Scanner scanner = new Scanner(System.in);

		Integer k = scanner.nextInt();
		Integer t = scanner.nextInt();
		Integer starts = scanner.nextInt();

		List<String> dna = new ArrayList<>();

		while(scanner.hasNext()){
			dna.add(scanner.next());
		}

		List<String> res = MotifEnumeration.gibbsSampler(dna,k,t,starts,100);

		for(String r : res){
			System.out.println(r);
		}

		System.out.println(MotifEnumeration.computeScore(res,k));

	}

}