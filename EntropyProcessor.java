import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf; //erfc function for Monobit test
import org.jtransforms.fft.DoubleFFT_1D; //Fourier Transform function for NIST test

public class EntropyProcessor {
	
	public static double ZALPHA = 2.5758293035489008;
	
	/*
	 * NOTE: Any additional information we need on the mouse inputs (timing for instance) can be added to the Pair class
	 */
	public static void main(String[] args) {
		if (args.length <= 1 || args.length >= 4) {
			System.err.println("ERROR: Must supply filenames and confidence value as a command line argument");
			System.exit(-1);
		}
		if (args.length == 2) { // For when you just give an input file and confidence value
			String filename = args[0];
			double confidence = Double.parseDouble(args[1]);
			
			// Calculate Z value based on provided confidence
			ZALPHA = calculateZ(confidence);
			
			ArrayList<Pair> data = readData(filename);
			
			process(data);
			
			double entropyEstimate = estimate(data);
			
			estimateMostCommon(data);
			
			estimateCollision(pairToBits(data));
			
			//String entropy = mix(processedData);
			
			//boolean testResult = test(entropy, confidence);
			
			System.out.println(frequencyNIST(pairToBits(data), confidence));
			System.out.println(runsNIST(pairToBits(data), confidence));
			System.out.println(fourierNIST(pairToBits(data), confidence));
		}
		else { // For when you give a input file, confidence value, and output file
			String inputFilename = args[0];
			String outputFilename = args[2];
			
			ArrayList<Pair> data = readData(inputFilename);
			
			try {
				FileWriter writer = new FileWriter(outputFilename);
				writer.write(pairToBits(data));
				writer.close();
			} catch (IOException e) {
				System.out.println("Something went screwy writing to file.");
				e.printStackTrace();
			}
			
		}
		
		// Print test result or something
	}
	
	/*
	 * Calculates log base 2 of a given integer as a double.
	 * Implemented by Joseph
	 */
	public static double log2(double x) {
		return Math.log(x) / Math.log(2.0);
	}
	
	/*
	 * Reads the filename and creates pairs of integers to represent the mouse data from that file
	 * Implemented by Joseph
	 */
	public static ArrayList<Pair> readData(String filename) {
		ArrayList<Pair> pairs = new ArrayList<Pair>();
		try {
			File coordinates = new File(filename);
			Scanner reader = new Scanner(coordinates);
			Pattern p = Pattern.compile("\\d+");
			Matcher m = null;
			int x = 0, y = 0;
			while (reader.hasNextLine()) {
				String line = reader.nextLine();
				m = p.matcher(line);
				m.find();
				x = Integer.parseInt(m.group());
				m.find();
				y = Integer.parseInt(m.group());
				pairs.add(new Pair(x, y));
			}
			reader.close();
		} catch (FileNotFoundException e) {
			System.err.println("ERROR: coordinates.txt not found");
			System.exit(-1);
		}
		return pairs;
	}
	
	/*
	 * Takes an array list of pairs of integers representing (X, Y) mouse coordinates and processes the data
	 * For instance, removing pairs at (0, 0).
	 * Returns a processed ArrayList of pairs of xintegers representing (X, Y) mouse coordinates
	 */
	public static ArrayList<Pair> process(ArrayList<Pair> data) {
		for (int i = 0; i < data.size(); i++) {
			if (data.get(i).x == 960 & data.get(i).y == 540) {
				data.remove(i);
			}
		}
		
		return data;
	}
	
	/*
	 * Taken from https://stats.stackexchange.com/questions/71788/percentile-to-z-score-in-php-or-java
	 */
	public static double calculateZ(double confidence) {
		double z = Math.sqrt(2) * Erf.erfcInv(2*confidence);
		System.out.println(z);
		return z;
	}
	
	/*
	 * Takes an ArrayList of pairs of integers representing (X, Y) mouse coordinates and returns a double representing the estimated amount of entropy from the data
	 * Implemented by Joseph
	 */
	public static double estimate(ArrayList<Pair> data) {
		// Count frequency of x, y pairs
		HashMap<Pair, Integer> frequency = new HashMap<Pair, Integer>();
		for (int i = 0; i < data.size(); i++) {
			if (frequency.containsKey(data.get(i))) {
				frequency.put(data.get(i), frequency.get(data.get(i)) + 1);
			}
			else {
				frequency.put(data.get(i), 1);
			}
		}
		// Sum the frequency of all pairs multiplied by their log(frequency) and multiply by negative one to calculate shannon entropy
		double sum = 0.0;
		double highestFrequency = 0.0;
		double probabilityOfPair = 0.0;
		
		// Take the amount of occurrences of each pair (i) and divide it by the total size of the data.  Then add it to the sum of shannon entropy by multiplying it by its log2
		for (Integer i : frequency.values()) {
			probabilityOfPair = (double)i / (double)data.size();
			if (probabilityOfPair > highestFrequency) {
				highestFrequency = probabilityOfPair;
			}
			sum += probabilityOfPair * log2(probabilityOfPair);
		}
		
		sum *= -1;
		System.out.println("Entropy Estimate = " + sum);
		
		double minEntropy = -1 * log2(highestFrequency);
		
		System.out.println("minEntropy = " + minEntropy);
		
		return sum;
	}
	
	/*
	 * 6.3.1 most common value estimate
	 * Computes the entropy of the sample based on the frequency of the most common sample and the P value
	 * Adapted code from https://github.com/usnistgov/SP800-90B_EntropyAssessment/blob/master/cpp/shared/most_common.h
	 * Implemented by Joseph
	 */
	public static double estimateMostCommon(ArrayList<Pair> data) {
		HashMap<Pair, Integer> frequency = new HashMap<Pair, Integer>();
		for (int i = 0; i < data.size(); i++) {
			if (frequency.containsKey(data.get(i))) {
				frequency.put(data.get(i), frequency.get(data.get(i)) + 1);
			}
			else {
				frequency.put(data.get(i), 1);
			}
		}
		
		int highestOccurences = 0;
		
		for (Integer i : frequency.values()) {
			if (i > highestOccurences) {
				highestOccurences = i;
			}
		}
		
		double highestFrequency = (double)highestOccurences / (double)data.size();
		
		double upperBound = Math.min(1.0, ZALPHA*Math.sqrt(highestFrequency*(1.0 - highestFrequency) / ((double)data.size() - 1)));
		
		double minEntropy = -log2(upperBound);
		System.out.println("Estimated min-entropy from estimateMostCommon = " + upperBound);
		
		return minEntropy;
	}
	
	/*
	 * 6.3.2 Collision Estimate
	 * Adapted code from https://github.com/usnistgov/SP800-90B_EntropyAssessment/blob/master/cpp/non_iid/collision_test.h
	 * Implemented by Joseph
	 */
	public static double estimateCollision(String binaryData) {
		int v, i ,j;
		int t_v;
		double X, s, p, lastP, pVal;
		double lvalue, hvalue;
		double hbound, lbound;
		double hdomain, ldomain;
		double entEst;
		
		i = 0;
		v = 0;
		s = 0.0;
		
		// compute wait times until collisions
		while(i < binaryData.length() - 1) {
			if (binaryData.charAt(i) == binaryData.charAt(i+1)) t_v = 2;
			else if (i < binaryData.length() - 2) t_v = 3;
			else break;
			
			v++;
			s += t_v*t_v;
			i += t_v;
		}
		
		// X is mean of t_v's, s is sample stdev, where
		// s^2 = (sum(t_v^2) - sum(t_v)^2/v) / (v-1)
		X = i / (double)v;
		s = Math.sqrt((s - (i*X)) / (v-1));
		
		X -= ZALPHA * s/Math.sqrt(v);
		
		// 2 is the smallest meaningful value here
		if (X < 2.0) X = 2.0;
		
		if (X < 2.5) {
			p = 0.5 + Math.sqrt(1.25 - 0.5 * X);
			entEst = -1 * log2(p);
			System.out.println("Collision estimate found");
		} else {
			System.out.println("Collision estimate could not find p, proceeding with lower bound for p");
			p = 0.5;
			entEst = 1.0;
		}
		
		System.out.println("Collision estimate = " + entEst);
		return entEst;
	}
	
	/*
	 * Takes an ArrayList of pairs of integers representing (X, Y) mouse coordinates and returns the string of entropy generated from this
	 * NOTE: String return type is a placeholder, return whatever data type would be most appropriate for representing entropy
	 */
	public static String mix(ArrayList<Pair> data) {
		int output = 0;
		int[] input = pairToBits2(data);
		int len = data.size();
		int shift = 0;

		for (int i=0; i < len -3; i++){
			int a = input[i];
			int b = input[i+1];

			input[i] = a + b;
			input[i+1] = Integer.rotateRight(b, shift);

			a = input[i];
			b = input[i+1];

			input[i+3] = a ^ b;

			shift = shift + 3;
			output = input[i+3];
		}
		return Integer.toBinaryString(output);
	}
	
	/*
	 * Takes a list of entropy strings and a confidence value and applies NIST tests to them to determine if the strings are random/unpredictable up to the confidence level
	 */
	public static boolean test(String entropy, double confidence) {
		frequencyNIST(entropy, confidence);
		runsNIST(entropy, confidence);
		fourierNIST(entropy, confidence);
		return false;
	}
	
	// Function to transform pairs of bits to a concatenated string of their binary equivalent.
	// Implemented by Aleks
	public static String pairToBits(ArrayList<Pair> data) {
		String bits = new String();
		
		int length = data.size();
		
		for (int i = 0; i < length; i++) {
			bits += Integer.toBinaryString(data.get(i).x);
			bits += Integer.toBinaryString(data.get(i).y);
		}
		return bits;
	}

	// Function which takes pair input and convert it to list of X+Y integers
	public static int[] pairToBits2(ArrayList<Pair> data) {
		String bits = new String();
		int length = data.size();
		int[] out = new int[length];

		for (int i = 0; i < length; i++) {
			bits = "";
			bits += Integer.toString(data.get(i).x);
			bits += Integer.toString(data.get(i).y);
			out[i] = Integer.parseInt(bits);
		}

		return out;
	}
	
	// Implementation of the NIST Frequency (Monobit) Test
	// Implemented by Aleks
	// Takes the bit string and a confidence/ significance value
	public static boolean frequencyNIST(String bits, double confidence) {
		System.out.println("NIST Frequency Test:");
		int length = bits.length();
		int sum = 0;
		
		int ones = 0;
		int zeroes = 0;
		
		for (int i = 0; i < length; i++) {
			sum += 2 * Character.getNumericValue(bits.charAt(i)) - 1;
			if (Character.getNumericValue(bits.charAt(i)) == 0) {
				zeroes++;
			}
			else {
				ones++;
			}
		}
		System.out.println("Sum: " + sum);
		System.out.println("Zeroes: " + zeroes);
		System.out.println("Ones: " + ones);
		
		double S = Math.abs(sum)/Math.sqrt(length);
		
		System.out.println("S: " + S);
		
		double pvalue = Erf.erfc(S/Math.sqrt(2));
		
		if (pvalue < confidence) {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
	public static boolean runsNIST(String bits, double confidence) {
		System.out.println("NIST Runs Test:");
		int length = bits.length();
		int V;
		int numOnes = 0;
		
		for (int i = 0; i < length; i++) {
			if (Character.getNumericValue(bits.charAt(i)) == 1) {
				numOnes++;
			}
		}
		
		double pi = numOnes/length;
		
		if (Math.abs(pi - 0.5) >= 2/Math.sqrt(length)) {
			System.out.println("Does not pass frequency test.");
			return false;
		}
		else {
			V = 1;
			for (int i = 1; i < length; i++) {
				if (Character.getNumericValue(bits.charAt(i)) != Character.getNumericValue(bits.charAt(i-1))) {
					V++;
				}
			}
		}
		
		double pvalue = Erf.erfc(Math.abs(V - 2 * length * pi * (1-pi)) / (2 * pi * (1-pi) * Math.sqrt(2*length)));
		
		if (pvalue < confidence) {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
	// Based on https://github.com/stamfest/randomtests/blob/master/src/main/java/net/stamfest/randomtests/nist/DiscreteFourierTransform.java
	public static boolean fourierNIST(String bits, double confidence) {
		System.out.println("NIST Discrete Fourier Transform Test:");
		int length = bits.length();
		double d;
		
		double X[] = new double[length];
		
		for (int i = 0; i < length; i++) {
            X[i] = 2 * (int) Character.getNumericValue(bits.charAt(i)) - 1;
        }
		
		DoubleFFT_1D fft = new DoubleFFT_1D(length);
		fft.realForward(X);
		
		double M[] = new double[length / 2 + 1]; // Not sure if that "+ 1" needs to be there.
		
		M[0] = Math.sqrt(X[0] * X[0]);
		M[length/2] = Math.sqrt(X[1] * X[1]);
		
		for (int i = 0; i < length / 2 - 1; i++) {
            M[i + 1] = Math.hypot(X[2 * i + 2], X[2 * i + 3]);
        }
		
		int count = 0;
        double upperBound = Math.sqrt(Math.log(1/0.05) * length);
        for (int i = 0; i < length / 2; i++) {
            if (M[i] < upperBound) {
                count++;
            }
        }
        
        double N1 = count;
        double N0 = 0.95 * length / 2.0;
        
        d = (N1 - N0) / Math.sqrt(length * 0.95 * 0.05 / 4);
        
        double pvalue = Erf.erfc(Math.abs(d) / Math.sqrt(2));
        
        if (pvalue < confidence) {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
	public static boolean cusumNIST (String bits, int mode, double confidence) {
		System.out.println("NIST Cumulative Sum Test:");
		int z;
        int zrev;
        int n = bits.length();
        double sqrtN = Math.sqrt(n);
        int S, sup, inf, k;
        double sum1, sum2;
        
        final NormalDistribution nd = new NormalDistribution();

        z = zrev = 0;
        S = 0;
        sup = 0;
        inf = 0;
        for (k = 0; k < n; k++) {
            if (Character.getNumericValue(bits.charAt(k)) != 0) {
                S++;
            } else {
                S--;
            }
            if (S > sup) {
                sup++;
            }
            if (S < inf) {
                inf--;
            }
            z = (sup > -inf) ? sup : -inf;
            zrev = (sup - S > S - inf) ? sup - S : S - inf;
        }
        
        double pvalue;

        if (mode == 0) {
	        // forward
	        sum1 = 0.0;
	        for (k = (-n / z + 1) / 4; k <= (n / z - 1) / 4; k++) {
	            sum1 += nd.cumulativeProbability(((4 * k + 1) * z) / sqrtN);
	            sum1 -= nd.cumulativeProbability(((4 * k - 1) * z) / sqrtN);
	        }
	        sum2 = 0.0;
	        for (k = (-n / z - 3) / 4; k <= (n / z - 1) / 4; k++) {
	            sum2 += nd.cumulativeProbability(((4 * k + 3) * z) / sqrtN);
	            sum2 -= nd.cumulativeProbability(((4 * k + 1) * z) / sqrtN);
	        }
	
	        pvalue = 1.0 - sum1 + sum2;
        }
        else {
	        // backwards
	        sum1 = 0.0;
	        for (k = (-n / zrev + 1) / 4; k <= (n / zrev - 1) / 4; k++) {
	            sum1 += nd.cumulativeProbability(((4 * k + 1) * zrev) / sqrtN);
	            sum1 -= nd.cumulativeProbability(((4 * k - 1) * zrev) / sqrtN);
	        }
	        sum2 = 0.0;
	        for (k = (-n / zrev - 3) / 4; k <= (n / zrev - 1) / 4; k++) {
	            sum2 += nd.cumulativeProbability(((4 * k + 3) * zrev) / sqrtN);
	            sum2 -= nd.cumulativeProbability(((4 * k + 1) * zrev) / sqrtN);
	        }
	        pvalue = 1.0 - sum1 + sum2;
        }
		
        if (pvalue < confidence) {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("P-value: " + pvalue);
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
}
