import java.io.*;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.math3.special.Erf; //erfc function for Monobit test, Gamma for block tests

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
			
			ArrayList<Integer> bitsArray = pairToBitsArray(data);
			String bitString			 = bitArrayToString(bitsArray);
			
			// For use with ENT program
			writeBytesToFile(bitString);
			
			double entropyEstimate 		 = estimate(bitsArray);
			
			estimateMostCommon(bitsArray);
			
			estimateCollision(bitArrayToString(bitsArray));
			
			//String entropy = mix(processedData);
			
			//boolean testResult = test(entropy, confidence);
			
			System.out.println("bitString.length() = " + bitString.length());
            
			System.out.println(NISTTests.frequencyNIST(bitArrayToString(pairToBitsArray(data)), confidence));
			System.out.println(NISTTests.blockFreqNIST(bitArrayToString(pairToBitsArray(data)), 3, confidence));
			System.out.println(NISTTests.runsNIST(bitArrayToString(pairToBitsArray(data)), confidence));
			System.out.println(NISTTests.blockRunsNIST(bitArrayToString(pairToBitsArray(data)), confidence));
			//System.out.println(); // Binary Matrix Rank Test
			System.out.println(NISTTests.fourierNIST(bitArrayToString(pairToBitsArray(data)), confidence));
			//System.out.println(); // Non-overlapping Template Matching Test
			//System.out.println(); // Overlapping Template Matching Test
			System.out.println(NISTTests.universalNIST(bitArrayToString(pairToBitsArray(data)), bitArrayToString(pairToBitsArray(data)).length()/5, 4, confidence));
			//System.out.println(); // Linear Complexity Test
			System.out.println(NISTTests.serialNIST(bitArrayToString(pairToBitsArray(data)), 3, confidence));
			System.out.println(NISTTests.approxEntropyNIST(bitArrayToString(pairToBitsArray(data)), 3, confidence));
			System.out.println(NISTTests.cusumNIST(bitArrayToString(pairToBitsArray(data)), 0, confidence));
			System.out.println(NISTTests.randomExcursionsNIST(bitArrayToString(pairToBitsArray(data)), confidence));
			System.out.println(NISTTests.randomExcursionsVariantNIST(bitArrayToString(pairToBitsArray(data)), confidence));
		}
	}
	
	/*
	 * Calculates log base 2 of a given integer as a double.
	 * Implemented by Joseph
	 */
	public static double log2(double x) {
		return Math.log(x) / Math.log(2.0);
	}
	
	
	/*
	 * Converts data to bytes for ENT program testing
	 * Implemented by Joseph
	 */
	public static void writeBytesToFile(String output) {
		BitSet bitSet = new BitSet(output.length());
		int counter = 0;
		for (char c : output.toCharArray()) {
			if (c == '1') {
				bitSet.set(counter);
			}
			counter++;
		}
		
		try (FileOutputStream stream = new FileOutputStream("output")) {
			stream.write(bitSet.toByteArray());
		} catch (Exception e) {
			System.out.println("ERROR: Error writing bytes to file.");
		}
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
		return z;
	}
	
	/*
	 * Takes an ArrayList of integers representing mouse coordinates and returns a double representing the estimated amount of entropy from the data
	 * Implemented by Joseph
	 */
	public static double estimate(ArrayList<Integer> data) {
		// Count frequency of integers from data
		HashMap<Integer, Integer> frequency = new HashMap<Integer, Integer>();
		
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
	public static double estimateMostCommon(ArrayList<Integer> bitsArray) {
		HashMap<Integer, Integer> frequency = new HashMap<Integer, Integer>();
		for (int i = 0; i < bitsArray.size(); i++) {
			if (frequency.containsKey(bitsArray.get(i))) {
				frequency.put(bitsArray.get(i), frequency.get(bitsArray.get(i)) + 1);
			}
			else {
				frequency.put(bitsArray.get(i), 1);
			}
		}
		
		int highestOccurences = 0;
		
		for (Integer i : frequency.values()) {
			if (i > highestOccurences) {
				highestOccurences = i;
			}
		}
		
		double highestFrequency = (double)highestOccurences / (double)bitsArray.size();
		
		double upperBound = Math.min(1.0, ZALPHA*Math.sqrt(highestFrequency*(1.0 - highestFrequency) / ((double)bitsArray.size() - 1)));
		
		double minEntropy = -log2(upperBound);
		System.out.println("Estimated min-entropy from estimateMostCommon = " + minEntropy);
		
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
		X = (double)i / (double)v;
		s = Math.sqrt((s - ((double)i*X)) / ((double)v-1));
		
		X -= ZALPHA * s/Math.sqrt((double)v);
		
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
		
		System.out.println("Collision Estimate: p = " + p);
		System.out.println("Collision Estimate: min entropy =  " + entEst);
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
		NISTTests.frequencyNIST(entropy, confidence);
		NISTTests.runsNIST(entropy, confidence);
		NISTTests.fourierNIST(entropy, confidence);
		return false;
	}
	
	/*
	 * Converts list of integers from pairToBitsArray into a string for use in the NIST tests
	 * Implemented by Joseph
	 */
	public static String bitArrayToString(ArrayList<Integer> bitarray) {
		String bitstring = "";
		for (int i = 0; i < bitarray.size(); i++) {
			bitstring += String.format("%10s", Integer.toBinaryString(bitarray.get(i)).replace(" ", "0"));
			//System.out.println(bitstring);
		}
		return bitstring;
	}
	
	/*
	 * Converts list of pairs to an arraylist of integers representing the 7 least significant bits of the X coordinate concatenated with the 3 least significant bits of the Y coordinate
	 * Implemented by Aleks & Joseph
	 */
	public static ArrayList<Integer> pairToBitsArray(ArrayList<Pair> data) {
		ArrayList<Integer> bits = new ArrayList<Integer>();
		
		int length = data.size();
		int xPrime = 0;
		int yPrime = 0;
		
		String xStrMask = "0000000001111111";
		int xMask = Integer.parseUnsignedInt(xStrMask, 2);
		
		String yStrMask = "0000000000000111";
		int yMask = Integer.parseUnsignedInt(yStrMask, 2);
		
		//System.out.println("xMask = " + xMask + "\t" + xStrMask);
		//System.out.println("yMask = " + yMask + "\t" + yStrMask);
		
		for (int i = 0; i < data.size(); i++) {
			xPrime = data.get(i).x;
			yPrime = data.get(i).y;
			
			//System.out.println(String.format("x:\t\t%16s", Integer.toBinaryString(xPrime)).replace(" ", "0"));
			//System.out.println(String.format("y:\t\t%16s", Integer.toBinaryString(yPrime)).replace(" ", "0"));
			
			xPrime &= xMask;
			yPrime &= yMask;
			
			String xString = String.format("%7s", Integer.toBinaryString(xPrime)).replace(" ", "0");
			String yString = String.format("%3s", Integer.toBinaryString(yPrime)).replace(" ", "0");
			//System.out.println(String.format("xPostMask:\t%16s", xString));
			//System.out.println(String.format("yPostMask:\t%16s", yString));
			//System.out.println(Integer.parseUnsignedInt(xString + yString, 2));
			bits.add(Integer.parseUnsignedInt(xString + yString, 2));
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
}
