import MessagesInDNA.FrequentWords;
import MessagesInDNA.FindClumps;
import java.util.Scanner;
import java.util.HashSet;

class Main{

	public static void main(String[] args){

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

	}

}