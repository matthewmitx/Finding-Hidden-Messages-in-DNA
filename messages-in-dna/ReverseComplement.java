package MessagesInDNA;
/* Program to find the reverse complement of a given DNA strand 
   Input: String segment
   Output: String reverseCompliment
*/

public class ReverseComplement{

	public static String reverseComplement(String segment){

		StringBuffer sb = new StringBuffer();

		for(char c : segment.toCharArray()){
			if(c == 'A')
				sb.append('T');
			else if(c == 'T')
				sb.append('A');
			else if(c == 'G')
				sb.append('C');
			else if(c == 'C')
				sb.append('G');
		}

		sb = sb.reverse();
		return sb.toString();
	}

}