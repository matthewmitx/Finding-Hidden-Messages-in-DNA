import MessagesInDNA.PatternMatching;

import java.util.*;

class Main{

	/* Main used for Approximate Pattern Matching program */
	public static void main(String[] args){

		Scanner scanner = new Scanner(System.in);
		String pattern = scanner.next();
		String segment = scanner.next();
		Integer d = scanner.nextInt();

		Integer count = PatternMatching.approxPatternCount(segment,pattern,d);

		System.out.println(count);

	}

	/* Main used for minimumSkew program */
	/* public static void main(String[] args){

		Scanner scanner = new Scanner(System.in);
		String genome = scanner.next();

		ArrayList<Integer> minimums = MinimumSkew.minimumSkew(genome);

		System.out.println(minimums.toString());

	} */

	/* Main used for FindClumps program */
	/* public static void main(String[] args){

		Scanner scanner = new Scanner(System.in);
		String segment = scanner.next();
		int k = scanner.nextInt();
		int l = scanner.nextInt();
		int t = scanner.nextInt();
		scanner.close();

		HashSet<String> clumps = FindClumps.findClumps(segment,k,l,t);

		for(String clump : clumps){
			System.out.println(clump);
		}

		System.out.println(clumps.size());

	} */

}