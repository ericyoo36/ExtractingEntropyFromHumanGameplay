import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;
import org.jtransforms.fft.DoubleFFT_1D;

// Implements 11 of the 15 NIST tests for testing randomness.

// Based on Java implementations by Peter Stamfest
// https://github.com/stamfest/randomtests/tree/master/src/main/java/net/stamfest/randomtests/nist

public class NISTTests {
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
	
	public static boolean blockFreqNIST (String bits, int blockLen, double confidence) {
		System.out.println("NIST Frequency Test within a Block:");
		double chi_squared;
        int blockCount;

        int n = bits.length();
        int discarded = n % blockLen;

        int i, j, blockSum;
        double sum, pi, v;

        blockCount = n / blockLen;
        /* # OF SUBSTRING BLOCKS      */
        sum = 0.0;

        for (i = 0; i < blockCount; i++) {
            blockSum = 0;
            int offset = i * blockLen;
            for (j = 0; j < blockLen; j++) {
                blockSum += Character.getNumericValue(bits.charAt((offset + j)));
            }
            pi = (double) blockSum / (double) blockLen;
            v = pi - 0.5;
            sum += v * v;
        }
        chi_squared = 4.0 * blockLen * sum;
        double pvalue = Gamma.regularizedGammaQ(blockCount/2.0, chi_squared/2.0);
       
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
	
	public static boolean blockRunsNIST (String bits, double confidence) {
		System.out.println("NIST Test for the Longest Run of Ones in a Block:");
		int blockLen;
        int numberOfBlocks;
        double chi2;
        int K;
//        boolean sequenceTooShort;

        int length = bits.length();
        double pi[] = new double[7];
        int run, v_n_obs, i, j, V[] = new int[7];

        int [] nu = new int[]{ 0, 0, 0, 0, 0, 0, 0 };

//        sequenceTooShort = false;
        if (length < 128) {
//            sequenceTooShort = true;
        	System.out.println("Sequence too short.");
            return false;
        }
        if (length < 6272) {
            K = 3;
            blockLen = 8;
            V[0] = 1;
            V[1] = 2;
            V[2] = 3;
            V[3] = 4;
            pi[0] = 0.21484375;
            pi[1] = 0.3671875;
            pi[2] = 0.23046875;
            pi[3] = 0.1875;
        } else if (length < 750000) {
            K = 5;
            blockLen = 128;
            V[0] = 4;
            V[1] = 5;
            V[2] = 6;
            V[3] = 7;
            V[4] = 8;
            V[5] = 9;
            pi[0] = 0.1174035788;
            pi[1] = 0.242955959;
            pi[2] = 0.249363483;
            pi[3] = 0.17517706;
            pi[4] = 0.102701071;
            pi[5] = 0.112398847;
        } else {
            K = 6;
            blockLen = 10000;
            V[0] = 10;
            V[1] = 11;
            V[2] = 12;
            V[3] = 13;
            V[4] = 14;
            V[5] = 15;
            V[6] = 16;
            pi[0] = 0.0882;
            pi[1] = 0.2092;
            pi[2] = 0.2483;
            pi[3] = 0.1933;
            pi[4] = 0.1208;
            pi[5] = 0.0675;
            pi[6] = 0.0727;
        }

        numberOfBlocks = length / blockLen;
        for (i = 0; i < numberOfBlocks; i++) {
            v_n_obs = 0;
            run = 0;
            for (j = 0; j < blockLen; j++) {
                if (Character.getNumericValue(bits.charAt((i * blockLen + j))) == 1) {
                    run++;
                    if (run > v_n_obs) {
                        v_n_obs = run;
                    }
                } else {
                    run = 0;
                }
            }
            if (v_n_obs < V[0]) {
                nu[0]++;
            }
            for (j = 0; j <= K; j++) {
                if (v_n_obs == V[j]) {
                    nu[j]++;
                }
            }
            if (v_n_obs > V[K]) {
                nu[K]++;
            }
        }

        chi2 = 0.0;
        for (i = 0; i <= K; i++) {
            chi2 += ((nu[i] - numberOfBlocks * pi[i]) * (nu[i] - numberOfBlocks * pi[i])) / (numberOfBlocks * pi[i]);
        }
        
        double pvalue = Gamma.regularizedGammaQ((double) (K / 2.0), chi2 / 2.0);
        
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

/*
	public static boolean rankNIST (String bits, double confidence) {
		
	}
	
	public static boolean templateMatchingNIST (String bits, int templateLength, double confidence) {
		
	}
*/

	public static boolean universalNIST (String bits, int numBlocks, int blockLength, double confidence) {
		System.out.println("NIST Universal Statistical Test:");
		int length = bits.length();

        int i, j, p;
        double arg, sqrt2, c;
        long T[];
        int decRep;
        
        final double variance[] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
                3.401, 3.410, 3.416, 3.419, 3.421 };

        final double expected_value[] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
                      8.1764248, 9.1723243, 10.170032, 11.168765,
                      12.168070, 13.167693, 14.167488, 15.167379 };

        int L = 5;
        if (length >= 387840) {
            L = 6;
        }
        if (length >= 904960) {
            L = 7;
        }
        if (length >= 2068480) {
            L = 8;
        }
        if (length >= 4654080) {
            L = 9;
        }
        if (length >= 10342400) {
            L = 10;
        }
        if (length >= 22753280) {
            L = 11;
        }
        if (length >= 49643520) {
            L = 12;
        }
        if (length >= 107560960) {
            L = 13;
        }
        if (length >= 231669760) {
            L = 14;
        }
        if (length >= 496435200) {
            L = 15;
        }
        if (length >= 1059061760) {
            L = 16;
        }

        int Q = 10 * (int) (1 << L);
        int K = (int) (Math.floor(length / L) - (double) Q);
        /* BLOCKS TO TEST */

        p = (int) 1 << L; // pow(2, L);

        T = new long[p];

        if ((L < 6) || (L > 16) || ((double) Q < 10 * (1 << L))) {
        	System.out.println("Poor conditions for test.");
            return false;
        }

        /* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
        c = 0.7 - 0.8 / (double) L + (4 + 32 / (double) L) * Math.pow(K, -3 / (double) L) / 15;
        double sigma = c * Math.sqrt(variance[L] / (double) K);
        sqrt2 = Math.sqrt(2);
        double log2 = Math.log(2);
        double sum = 0.0;

        for (i = 0; i < p; i++) {
            T[i] = 0;
        }
        for (i = 1; i <= Q; i++) {
            /* INITIALIZE TABLE */
            decRep = 0;
            for (j = 0; j < L; j++) {
                decRep += Character.getNumericValue(bits.charAt((i - 1) * L + j)) * (long) Math.pow(2, L - 1 - j);
            }
            T[decRep] = (long) i;
        }
        for (i = Q + 1; i <= Q + K; i++) {
            /* PROCESS BLOCKS */
            decRep = 0;
            for (j = 0; j < L; j++) {
                decRep += Character.getNumericValue(bits.charAt((i - 1) * L + j)) * (long) Math.pow(2, L - 1 - j);
            }
            sum += Math.log((double) i - T[decRep]) / log2;
            T[decRep] = i;
        }
        double phi = (double) (sum / (double) K);
        arg = Math.abs(phi - expected_value[L]) / (sqrt2 * sigma);
        
        double pvalue = Erf.erfc(arg);
        
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
	
	public static boolean serialNIST (String bits, int blockLength, double confidence) {
		double psim0, psim1, psim2, del1, del2;
        int length = bits.length();

        psim0 = psi2(blockLength, bits);
        psim1 = psi2(blockLength - 1, bits);
        psim2 = psi2(blockLength - 2, bits);
        del1 = psim0 - psim1;
        del2 = psim0 - 2.0 * psim1 + psim2;
        double pvalue1 = Gamma.regularizedGammaQ(Math.pow(2, blockLength - 1) / 2, del1 / 2.0);
        double pvalue2 = Gamma.regularizedGammaQ(Math.pow(2, blockLength - 2) / 2, del2 / 2.0);
	
        if (pvalue1 < confidence && pvalue2 < confidence) {
			System.out.println("P-value 1: " + pvalue1);
			System.out.println("P-value 2: " + pvalue2);
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("P-value 1: " + pvalue1);
			System.out.println("P-value 2: " + pvalue2);
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
	// psi2 function for serialNIST
	private static double psi2(int m, String bits) {
        int length = bits.length();

        int i, j, k, powLen;
        double sum, numOfBlocks;
        int P[];

        if ((m == 0) || (m == -1)) {
            return 0.0;
        }
        numOfBlocks = length;
        powLen = (int) Math.pow(2, m + 1) - 1;

        P = new int[powLen];
        for (i = 1; i < powLen - 1; i++) {
            P[i] = 0;
            /* INITIALIZE NODES */
        }
        for (i = 0; i < numOfBlocks; i++) {
            /* COMPUTE FREQUENCY */
            k = 1;
            for (j = 0; j < m; j++) {
                if (Character.getNumericValue(bits.charAt(((i + j) % length))) == 0) {
                    k *= 2;
                } else if (Character.getNumericValue(bits.charAt(((i + j) % length))) == 1) {
                    k = 2 * k + 1;
                }
            }
            P[k - 1]++;
        }
        sum = 0.0;
        for (i = (int) Math.pow(2, m) - 1; i < (int) Math.pow(2, m + 1) - 1; i++) {
            sum += Math.pow(P[i], 2);
        }
        sum = (sum * Math.pow(2, m) / (double) length) - (double) length;

        return sum;

    }
	
	public static boolean approxEntropyNIST(String bits, int blockLength, double confidence) {
		int seqLength;
        double chi_squared;
        double ApEn[] = new double[2], apen;
       
        int i, j, k, r, blockSize, powLen, index;
        double sum, numOfBlocks;
        int P[];

        seqLength = bits.length();
        r = 0;

        for (blockSize = blockLength; blockSize <= blockLength + 1; blockSize++) {
            if (blockSize == 0) {
                ApEn[0] = 0.00;
                r++;
            } else {
                numOfBlocks = (double) seqLength;
                powLen = (int) Math.pow(2, blockSize + 1) - 1;
                P = new int[powLen];

                for (i = 1; i < powLen - 1; i++) {
                    P[i] = 0;
                }
                for (i = 0; i < numOfBlocks; i++) {
                    /* COMPUTE FREQUENCY */
                    k = 1;
                    for (j = 0; j < blockSize; j++) {
                        k <<= 1;
                        if ((int) Character.getNumericValue(bits.charAt((i + j) % seqLength)) == 1) {
                            k++;
                        }
                    }
                    P[k - 1]++;
                }
                /* DISPLAY FREQUENCY */
                sum = 0.0;
                index = (int) Math.pow(2, blockSize) - 1;
                for (i = 0; i < (int) Math.pow(2, blockSize); i++) {
                    if (P[index] > 0) {
                        sum += P[index] * Math.log(P[index] / numOfBlocks);
                    }
                    index++;
                }
                sum /= numOfBlocks;
                ApEn[r] = sum;
                r++;
            }
        }
        apen = ApEn[0] - ApEn[1];

        chi_squared = 2.0 * seqLength * (Math.log(2) - apen);
        double pvalue = Gamma.regularizedGammaQ(Math.pow(2, blockLength - 1), chi_squared / 2.0);
        
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
	
	public static boolean randomExcursionsNIST(String bits, double confidence) {
		int n = bits.length();
        int J;
        double constraint;

        int b, i, j, k, x;
        int cycleStart, cycleStop, cycle[], S_k[];
        int counter[] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        double sum, nu[][] = new double[6][8];

        S_k = new int[n];
        cycle = new int[Math.max(1000, n / 100)];

        J = 0;
        
        int stateX[] = { -4, -3, -2, -1, 1, 2, 3, 4 };
        double pi[][] = {
            { 0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000 },
            { 0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000 },
            { 0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625 },
            { 0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143 },
            { 0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051 } };
        
        /* DETERMINE CYCLES */
        S_k[0] = 2 * (int) Character.getNumericValue(bits.charAt(0)) - 1;
        for (i = 1; i < n; i++) {
            S_k[i] = S_k[i - 1] + 2 * Character.getNumericValue(bits.charAt(i)) - 1;
            if (S_k[i] == 0) {
                J++;
                if (J > Math.max(1000, n / 100)) {
                    throw new IndexOutOfBoundsException("ERROR IN randomExcursions: EXCEEDING THE MAX NUMBER OF CYCLES EXPECTED");
                }
                cycle[J] = i;
            }
        }
        
        if (S_k[n - 1] != 0) {
            J++;
        }
        
        cycle[J] = n;

        constraint = Math.max(0.005 * Math.pow(n, 0.5), 500);
        if (J < constraint) {
            /* The NIST publications section 3.14 indicates that tests should 
            be rejected if this constraint is not met. So we do this... */
           
        	System.out.println("Data does not fit constraint.");
        	return false;
        }
        

        cycleStart = 0;
        cycleStop = cycle[1];
        for (k = 0; k < 6; k++) {
            for (i = 0; i < 8; i++) {
                nu[k][i] = 0.;
            }
        }
        for (j = 1; j <= J; j++) {
            /* FOR EACH CYCLE */
            for (i = 0; i < 8; i++) {
                counter[i] = 0;
            }
            for (i = cycleStart; i < cycleStop; i++) {
                if ((S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1)) {
                    if (S_k[i] < 0) {
                        b = 4;
                    } else {
                        b = 3;
                    }
                    counter[S_k[i] + b]++;
                }
            }
            cycleStart = cycle[j] + 1;
            if (j < J) {
                cycleStop = cycle[j + 1];
            }

            for (i = 0; i < 8; i++) {
                if ((counter[i] >= 0) && (counter[i] <= 4)) {
                    nu[counter[i]][i]++;
                } else if (counter[i] >= 5) {
                    nu[5][i]++;
                }
            }
        }
        
        double pvalues[] = new double [8];

        for (i = 0; i < 8; i++) {
            x = stateX[i];
            sum = 0.;
            for (k = 0; k < 6; k++) {
                sum += Math.pow(nu[k][i] - J * pi[(int) Math.abs(x)][k], 2) / (J * pi[(int) Math.abs(x)][k]);
            }
            
            pvalues[i] = Gamma.regularizedGammaQ(2.5, sum / 2.0);
        }
        
        double pvalue;
        boolean pass = true;
        
        for (int iter = 0; iter < 8; iter++) {
        	pvalue = pvalues[iter];
        	
        	System.out.println("P-value [" + iter + "] : " + pvalue);
        	
	        if (pvalue < confidence) {
				pass = false;
			}
        }
        
        if (!pass) {
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
	public static boolean randomExcursionsVariantNIST(String bits, double confidence) {
		int length;
        int J;
        int constraint;
    
        length = bits.length();
        int i, p, x, count, S_k[];

        S_k = new int[length];

        J = 0;
        
        int stateX[] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
        
        S_k[0] = 2 * (int) Character.getNumericValue(bits.charAt(0)) - 1;
        for (i = 1; i < length; i++) {
            S_k[i] = S_k[i - 1] + 2 * Character.getNumericValue(bits.charAt(i)) - 1;
            if (S_k[i] == 0) {
                J++;
            }
        }
        if (S_k[length - 1] != 0) {
            J++;
        }
        
        double pvalues[] = new double[16];

        constraint = (int) Math.max(0.005 * Math.pow(length, 0.5), 500);
        if (J < constraint) {
        	return false;

        } 
        else {
            for (p = 0; p <= 17; p++) {
                x = stateX[p];
                count = 0;
                for (i = 0; i < length; i++) {
                    if (S_k[i] == x) {
                        count++;
                    }
                }
                
                pvalues[p] = Erf.erfc(Math.abs(count - J) / (Math.sqrt(2.0 * J * (4.0 * Math.abs(x) - 2))));
            }
        }
        
        double pvalue;
        boolean pass = true;
        
        for (int iter = 0; iter < 8; iter++) {
        	pvalue = pvalues[iter];
        	
        	System.out.println("P-value [" + iter + "] : " + pvalue);
        	
	        if (pvalue < confidence) {
				pass = false;
			}
        }
        
        if (!pass) {
			System.out.println("Null hypothesis rejected (Not random).");
			return false;
		}
		else {
			System.out.println("Null hypothesis accepted (Sufficiently random).");
			return true;
		}
	}
	
}
