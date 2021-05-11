package MessagesInDNA;

/* Helper functions for MessagesInDNA programs */

import java.util.*;
import java.lang.Math.*;

public class Helper{

	/* 	Generates and returns a random number 1 through n */
	public static Integer randomNumber(Integer n){
		double rand = Math.random();
		int wholeNum = (int) (rand * (double) n);
		return wholeNum;
	}

	/* 	Generates and returns a random number 1 through size of probabilities list 
	`	based on the the probability distribution contained in the list */
	public static Integer random(List<Double> probabilities){
		
		double sum = 0;
		for(Double n : probabilities){
			sum += n;
		}

		double rand = Math.random() * sum;
		double runningSum = probabilities.get(0);

		if(rand < runningSum)
			return 0;


		for(int i = 1; i < probabilities.size(); i++){
			runningSum += probabilities.get(i);
			if(rand < runningSum)
				return i;
		}
		return probabilities.size()-1;
	}


}