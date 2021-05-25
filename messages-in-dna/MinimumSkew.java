package MessagesInDNA;

/* Contains functions
		minimumSkew
*/

import java.util.ArrayList;

public class MinimumSkew{

	/* 	Function to find the possible location of the replication origin by
   		measuring the total difference between guanine and cytosine level at
   		each nucleotide.
   		Input: String genome
   		Output: ArrayList of integers that minimize the skew[i](genome)

   		The input string genome begins anywhere along the genome and is
   		directed in the 5' -> 3' direction.
	*/
	public static ArrayList<Integer> minimumSkew(String genome){

		int skew = 0;
		int minSkew = 0;

		ArrayList<Integer> minimums = new ArrayList<>();

		for(int i = 0; i < genome.length(); i++){
			if(genome.charAt(i) == 'G')
				skew++;
			else if(genome.charAt(i) == 'C')
				skew--;

			if(skew == minSkew){
				minimums.add(i);
			}
			else if(skew < minSkew){
				minimums.clear();
				minimums.add(i);
				minSkew = skew;
			}
		
		}
		return minimums;
	}

}