import java.util.ArrayList;

public class EntropyProcessor {
	
	/*
	 * NOTE: Any additional information we need on the mouse inputs (timing for instance) can be added to the Pair class
	 */

	public static void main(String[] args) {
		if (args.length != 3) {
			System.err.println("ERROR: Must supply filename and confidence value as a command line argument");
		}
		String filename = args[1];
		double confidence = Double.parseDouble(args[2]);
		
		ArrayList<Pair> data = readData(filename);
		
		ArrayList<Pair> processedData = process(data);
		
		double entropyEstimate = estimate(processedData);
		
		String[] entropy = mix(processedData);
		
		boolean testResult = test(entropy, confidence);
		
		// Print test result or something
	}
	
	/*
	 * Reads the filename and creates pairs of integers to represent the mouse data from that file
	 */
	public static ArrayList<Pair> readData(String filename) {
		return null;
	}
	
	/*
	 * Takes an array list of pairs of integers representing (X, Y) mouse coordinates and processes the data
	 * For instance, removing pairs at (0, 0).
	 * Returns a processed ArrayList of pairs of integers representing (X, Y) mouse coordinates
	 */
	public static ArrayList<Pair> process(ArrayList<Pair> data) {
		return data;
	}
	
	/*
	 * Takes an ArrayList of pairs of integers representing (X, Y) mouse coordinates and returns a double representing the estimated amount of entropy from the data
	 */
	public static double estimate(ArrayList<Pair> data) {
		return 0.0;
	}
	
	/*
	 * Takes an ArrayList of pairs of integers representing (X, Y) mouse coordinates and returns the string of entropy generated from this
	 * NOTE: String return type is a placeholder, return whatever data type would be most appropriate for representing entropy
	 */
	public static String[] mix(ArrayList<Pair> data) {
		String[] strings = new String[1];
		return strings;
	}
	
	/*
	 * Takes a list of entropy strings and a confidence value and applies NIST tests to them to determine if the strings are random/unpredictable up to the confidence level
	 */
	public static boolean test(String[] entropy, double confidence) {
		return false;
	}
}
