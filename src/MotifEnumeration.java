package MessagesInDNA;

/* Contains functions:
	motifEnumeration
	medianString
	distanceBetweenPatternAndStrings
	profileMostProbableKmer
	greedyMotifSearch
	createProfileMatrix
	createLaplaceProfileMatrix
	computeScore
	computeProbability
*/

import java.util.*;

public class MotifEnumeration{

	public static final Integer NUM_NUCLEOTIDES = 4;

	/* 	Function to find the all k-mers of specified size k from a given list of 
		DNA segments. The k-mers may contain at most d mismatches but they should
		be present within every DNA segment in list dna. 
  	   	Input:  String dna 
   			    Integer k 
   			    Integer d
  	   	Output: Set of all k-mers commonly found in each String in dna list
	*/
	public static Set<String> motifEnumeration(List<String> dna, Integer k, Integer d){

		List<Set<String>> patterns = new ArrayList<>();
		for(int i = 0; i < dna.size(); i++){
			patterns.add(new HashSet<String>());
		}

		for(int i = 0; i < dna.size(); i++){
			String upstreamRegion = dna.get(i);
			for(int j = 0; j < upstreamRegion.length()-k+1; j++){
				String pattern = upstreamRegion.substring(j,j+k);
				HashSet<String> neighbors = FrequentWords.findNeighbors(pattern,d);
				patterns.get(i).addAll(neighbors);
			}
		}

		Set<String> motifs = patterns.get(0);

		for(int i = 1; i < patterns.size(); i++){
			motifs.retainAll(patterns.get(i));
		}

		return motifs;
	}

	/* 	Function to find the single string Pattern that has the minimum distance
		from the list of dna segments - dna. It uses the auxillary function
		distanceBetweenPatternAndStrings to calculate the distance between each
		potential pattern and the input dna list.
  	   	Input:  Integer k (length of pattern)
   			    List<String> dna
  	   	Output: String pattern
	*/
	public static String medianString(Integer k, List<String> dna){

		StringBuilder generic = new StringBuilder();
		for(int i = 0; i < k; i++){
			generic.append("A");
		}

		Set<String> allPatterns = FrequentWords.findNeighbors(generic.toString(),k);

		Integer finalDistance = Integer.MAX_VALUE;
		String finalPattern = "";

		for(String currentPattern : allPatterns){
			int currentDistance = distanceBetweenPatternAndStrings(currentPattern,dna);
			if(currentDistance < finalDistance){
				finalDistance = currentDistance;
				finalPattern = currentPattern;
			}
		}

		return finalPattern;
	}


	/* 	Function to calcuate the distance between a given pattern and a list of dna
		strings. The distance between each dna and the pattern is equal to the minimum
		difference of characters between pattern and any substring in the dna of the 
		same length.
  	   	Input:  String pattern 
   			    List<String> dna (list of separate dna segments)
  	   	Output: Sum of differences between pattern and each segment in dna
	*/
	public static Integer distanceBetweenPatternAndStrings(String pattern, List<String> dna){

		int k = pattern.length();
		int totalDistance = 0;

		for(String dnaSegment : dna){
			int minimumDistance = k+1;
			for(int i = 0; i < dnaSegment.length()-k+1; i++){
				String currentPattern = dnaSegment.substring(i,i+k);
				int currentDistance = PatternMatching.hammingDistance(pattern,currentPattern);
				if(currentDistance < minimumDistance){
					minimumDistance = currentDistance;
				}
			}
			totalDistance += minimumDistance;
		}

		return totalDistance;
	}


	/* 
		Function to find the substring on size k from text that is most likely to 
		occur. Each substring's probability is measured by using the 4xk profile
		matrix that expresses the likelihood of each nucleotides occurrence at the
		k different positions.
		Input: 	String text
				Integer k
				float[][] profile (matrix of size 4xk)
		Output:	String most probabily substring contained in text of size k
	*/
	public static String profileMostProbableKmer(String text, Integer k, float[][] profile){

		String finalPattern = text.substring(0,k);
		float finalProbability = computeProbability(finalPattern,profile);

		for(int i = 1; i < text.length()-k+1; i++){

			String pattern = text.substring(i,i+k);
			float probability = computeProbability(pattern,profile);

			if(probability > finalProbability){
				finalProbability = probability;
				finalPattern = pattern;
			}
		}
		return finalPattern;
	}


	/* 
		Function to search for the best set of motifs from dna stands - dna. The algorithm 
		does this by choosing a set of motifs from all strands (2 to t) for each k-mer 
		existing in the first strand. The function then returns the best set of motifs.
		Input: 	List<String> dna
				Integer k
				Integer t (size of dna list)
		Output:	List of motifs (one from each dna strand) that are most probable

		Function can be optimized by eliminating the conversion between a list of strings
		and a 2d float array for the profile matricies.
	*/
	public static List<String> greedyMotifSearch(List<String> dna, Integer k, Integer t){

		List<String> bestMotifs = new ArrayList<>();

		for(String strand : dna){
			bestMotifs.add(strand.substring(0,k));
		}

		String baseStrand = dna.get(0);
		List<String> otherStrands = new ArrayList<>(dna);
		otherStrands.remove(0);

		for(int i = 0; i < baseStrand.length()-k+1; i++){
			// base motif from first strand to guide the choice of motifs from other strands
			String guideMotif = baseStrand.substring(i,i+k);
			List<String> allCurrMotifs = new ArrayList<>();
			allCurrMotifs.add(guideMotif);

			for(String strand : otherStrands){
				float[][] profileMatrix = createLaplaceProfileMatrix(allCurrMotifs,k);
				String nextMotif = profileMostProbableKmer(strand,k,profileMatrix);
				allCurrMotifs.add(nextMotif);
			}

			if(computeScore(createProfileMatrix(allCurrMotifs,k),k) < computeScore(createProfileMatrix(bestMotifs,k),k)){
				bestMotifs = allCurrMotifs;
			}
		}
		return bestMotifs;
	}


	/*
		Function to find a set of motifs from a set of strings dna. The motifs found will be of length k. The function
		uses a ranadomized approach in which it chooses a random pattern from each string in dna. It motifies this 
		pattern on each iteration by using the set of patterns' profile matrix to choose the patterns that better match.
		This algorithm works by running many times. Inevitably, the algorithm will choose the actual motifs on one of
		these iterations, therefore having the lowest score and giving the best set of motifs.
		Input: 	List<String> dna
				Integer k
				Integer iteration (number of times to run the algorithm)
		Output:	List<String> bestMotifs resulting from function run -iteration- times
	*/
	public static List<String> randomizedMotifSearch(List<String> dna, Integer k, Integer iterations){

		if(dna.size() <= 0)
			return new ArrayList<String>();

		int length = dna.get(0).length();

		List<String> absoluteBestMotifs = new ArrayList<>();
		int absoulteBestScore = Integer.MAX_VALUE;

		for(int i = 0; i < iterations; i++){

			List<String> motifs = new ArrayList<>();

			for(String strand : dna){
				int rand = Helper.randomNumber(length-k+1);
				String motif = strand.substring(rand,rand+k);
				motifs.add(motif);
			}

			List<String> bestMotifs = new ArrayList<>(motifs);

			while(true){
				float[][] profileMatrix = createLaplaceProfileMatrix(motifs,k);
				motifs = generateMostProbableMotifs(dna,k,profileMatrix);
				int bestScore = computeScore(bestMotifs,k);
				if(computeScore(motifs,k) < bestScore){
					bestMotifs = new ArrayList<>(motifs);
				}
				else{
					if(bestScore < absoulteBestScore){
						absoulteBestScore = bestScore;
						absoluteBestMotifs = new ArrayList<>(bestMotifs);
					}
					break;
				}
			}
		}
		return absoluteBestMotifs;
	}


	/*
		Function to determine the most probable candidates to be motifs in a set of strings dna, given
		a profile matrix. 
		Input: 	List<String> dna
				Integer k
				float[][] profile
		Output:	List<String> with the best candidates to be the motif strings
	*/
	public static List<String> generateMostProbableMotifs(List<String> dna, Integer k, float[][] profile){

		List<String> motifs = new ArrayList<>();

		for(String strand : dna){

			String bestMotif = strand.substring(0,k);
			double bestProbability = computeProbability(bestMotif,profile);

			double probability;

			for(int i = 1; i < strand.length()-k+1; i++){
				String motif = strand.substring(i,i+k);
				probability = computeProbability(motif,profile);

				if(probability > bestProbability){
					bestProbability = probability;
					bestMotif = motif;
				}
			}
			motifs.add(bestMotif);
		}
		return motifs;
	}


	/* 
		Function to create a profile matrix for a given set of motifs. This profile matrix
		contains the likelihood of each of the four nucleotides at each position from the
		set of motifs.
		Input: 	List<String> motifs
				Integer k
		Output:	4 x k matrix with probabilities of each nucleotide A,C,G,T
	*/
	public static float[][] createProfileMatrix(List<String> motifs, Integer k){

		if(motifs.isEmpty())
			return new float[0][0];

		float numMotifs = motifs.size();

		float[][] profileMatrix = new float[NUM_NUCLEOTIDES][k];

		for(int i = 0; i < k; i++){
			int a = 0,c = 0,g = 0,t = 0;
			for(int j = 0; j < numMotifs; j++){
				char curr = motifs.get(j).charAt(i);
				if(curr == 'A'){
					a++;
				}
				else if(curr == 'C'){
					c++;
				}
				else if(curr == 'G'){
					g++;
				}
				else if(curr == 'T'){
					t++;
				}
			}
			profileMatrix[0][i] = a / numMotifs;
			profileMatrix[1][i] = c / numMotifs;
			profileMatrix[2][i] = g / numMotifs;
			profileMatrix[3][i] = t / numMotifs;
		}

		return profileMatrix;
	}


	/* 
		Function to create a profile matrix for a given set of motifs based on the Laplace
		rule of succession. This profile matrix contains the likelihood of each of the four
		nucleotides at each position from the set of motifs. The Laplace rule ensure that no
		single nucleotide has a probability of zero at any position in the profile matrix by
		using pseduocounts for nucleotides that do not appear in the sample size.
		Input: 	List<String> motifs
				Integer k
		Output:	4 x k matrix with probabilities of each nucleotide A,C,G,T
	*/
	public static float[][] createLaplaceProfileMatrix(List<String> motifs, Integer k){

		if(motifs.isEmpty())
			return new float[0][0];

		float numMotifs = motifs.size();

		float[][] profileMatrix = new float[NUM_NUCLEOTIDES][k];

		for(int i = 0; i < k; i++){
			int a = 1,c = 1,g = 1,t = 1;
			for(int j = 0; j < numMotifs; j++){
				char curr = motifs.get(j).charAt(i);
				if(curr == 'A'){
					a++;
				}
				else if(curr == 'C'){
					c++;
				}
				else if(curr == 'G'){
					g++;
				}
				else if(curr == 'T'){
					t++;
				}
			}
			profileMatrix[0][i] = a / (numMotifs+4);
			profileMatrix[1][i] = c / (numMotifs+4);
			profileMatrix[2][i] = g / (numMotifs+4);
			profileMatrix[3][i] = t / (numMotifs+4);
		}

		return profileMatrix;
	}


	/* 
		Function to compute the score of a profile matrix. Score is similar to entropy 
		in that it attempts to measure the certainty of a specific outcome. 
		Input: 	float[][] profile
				Integer k
		Output:	Float score of the profile matrix
	*/
	public static Float computeScore(float[][] profile, Integer k){

		float sum = 0; 

		for(int i = 0; i < k; i++){
			float greatest = 0;
			float currSum = 0;
			for(int j = 0; j < NUM_NUCLEOTIDES; j++){
				currSum += profile[j][i];
				if(profile[j][i] > greatest)
					greatest = profile[j][i];
			}
			currSum -= greatest;
			sum += currSum;
		}
		return sum;
	}

	/* 	Same compute score function but uses a list of motifs as input rather than
		a profile matrix.
	*/
	public static Integer computeScore(List<String> motifs, Integer k){

		int sum = 0; 
		int[] counts = new int[NUM_NUCLEOTIDES];

		for(int i = 0; i < k; i++){
			for(int j = 0; j < motifs.size(); j++){
				switch(motifs.get(j).charAt(i)){
					case 'A' :
						counts[0]++;
						break;
					case 'C' :
						counts[1]++;
						break;
					case 'G' :
						counts[2]++;
						break;
					case 'T' :
						counts[3]++;
						break;
				}
			}
			int greatest = greatest(counts);
			sum += (motifs.size() - greatest);
			counts[0] = counts[1] = counts[2] = counts[3] = 0;
		}
		return sum;
	}

	/* Returns the greatest number in an array */
	public static Integer greatest(int[] nums){
		if(nums.length == 0)
			return 0;

		int greatest = nums[0];
		for(int i = 1; i < nums.length; i++){
			if(nums[i] > greatest){
				greatest = nums[i];
			}
		}
		return greatest;
	}


	/* 
		Function to calculate the probability of a string pattern's occurence 
		given a probabilty matrix profile.
		Input: 	String pattern
				float[][] profile (matrix of size 4xk)
		Output:	Integer probabilty of string's occurence
	*/
	public static Float computeProbability(String pattern, float[][] profile){

		float probability = 1;
		for(int i = 0; i < pattern.length(); i++){
			switch(pattern.charAt(i)){
				case 'A' :
					probability *= profile[0][i];
					break;
				case 'C' :
					probability *= profile[1][i];
					break;
				case 'G' :
					probability *= profile[2][i];
					break;
				case 'T' :
					probability *= profile[3][i];
					break;
			}
			if(probability == 0)
				break;
		}
		return probability;
	}



}