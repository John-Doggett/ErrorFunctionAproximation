import java.math.BigDecimal;
import java.math.MathContext;
import java.util.ArrayList;

public class Aproximate {
static private MathContext PRECISION;
	public static double[] aproximate(Function func, double lowerBoundd, double upperBoundd, int numPoints, int polynomialPower, int precision) {
		PRECISION = new MathContext(precision);
		BigDecimal lowerBound = new BigDecimal(lowerBoundd, PRECISION);
		BigDecimal upperBound = new BigDecimal(upperBoundd, PRECISION);

		BigDecimal increment = (upperBound.subtract(lowerBound, PRECISION)).divide(new BigDecimal(numPoints - 1, PRECISION), PRECISION); //((double) (numPoints - 1));
		BigDecimal[][] A = new BigDecimal[numPoints][polynomialPower];
		for (int a = 0; a < A.length; a++) {
			A[a][0] = new BigDecimal(1, PRECISION);
		}
		for (BigDecimal a = new BigDecimal(0, PRECISION); a.doubleValue() < A.length; a = a.add(new BigDecimal(1, PRECISION), PRECISION)) {
			A[(int) a.doubleValue()][1] = lowerBound.add( a.multiply(increment, PRECISION), PRECISION);
		}
		for (int a = 0; a < A.length; a++) {
			for (int b = 2; b < A[0].length; b++) {
				A[a][b] = A[a][1].pow(b);
			}
		}
		// A setup

		BigDecimal[][] B = new BigDecimal[A.length][1];
		for (BigDecimal a = new BigDecimal(0, PRECISION); a.doubleValue() < B.length; a = a.add(new BigDecimal(1, PRECISION), PRECISION)) {
			BigDecimal temp = func.evaulate(a.multiply(increment, PRECISION).add(lowerBound,PRECISION), precision);
			B[(int) a.doubleValue()][0] = temp;
		}
		// B setup

		BigDecimal[][] AT = transpose(A);
		BigDecimal[][] C = matrixMultiply(AT, A);
		BigDecimal[][] D = matrixMultiply(AT, B);
		BigDecimal[][] E = concat(C, D);
		BigDecimal[][] FINAL = RREF(E);
		double[] output = new double[FINAL.length];
		for(int a = 0; a < FINAL.length; a++) {
			double temp = FINAL[a][FINAL[a].length-1].doubleValue();
			output[a] = temp;
		}
		return output;
	}

	private static BigDecimal[][] RREF(BigDecimal[][] inputMatrix) {
		BigDecimal[][] input = new BigDecimal[inputMatrix.length][inputMatrix[0].length];
		for (int a = 0; a < inputMatrix.length; a++) {
			for (int b = 0; b < inputMatrix[0].length; b++) {
				input[a][b] = inputMatrix[a][b];
			}
		}

		ArrayList<Integer> pivotPositions = new ArrayList<>();
		int length = input[0].length;
		for (int a = 0; a < input.length; a++) {
			if (input[a].length != length) {
				throw new IllegalArgumentException("Matrix is not mxn!");
			}
		}

		// REF Forward Operation
		int rowAdjustor = 0;
		for (int a = 0; a < input[0].length && a - rowAdjustor < input.length; a++) {

			int largest = a - rowAdjustor;
			for (int b = a - rowAdjustor; b < input.length; b++) {
				//if (Math.abs(input[b][a]) > Math.abs(input[largest][a])) {
				if ((input[b][a]).abs(PRECISION).compareTo((input[largest][a]).abs(PRECISION)) > 0) {

					largest = b;
				}
			}
			swapRow(input, a - rowAdjustor, largest);
			if ((input[a - rowAdjustor][a]).doubleValue() == 0.0) {
				rowAdjustor++;
			} else {
				pivotPositions.add(a);

				for (int b = a + 1 - rowAdjustor; b < input.length; b++) {
					multiplyAddRow(input, a - rowAdjustor, b,
							new BigDecimal(-1.0, PRECISION).multiply(input[b][a].divide(input[a - rowAdjustor][a], PRECISION), PRECISION));
				}

			}
		}

		// BACKWARD SUBSTITUTION
		// We know where the pivot positions are, which we can use to ignore critical
		// double error.
		// We are in REF.
		for (int a = pivotPositions.size(); a < input.length; a++) {
			for (int b = 0; b < input[a].length; b++) {
				if (a != pivotPositions.size() - 1 && b != input[0].length - 1
						&& pivotPositions.get(pivotPositions.size() - 1) == input[0].length - 1) {
					input[a][b] = new BigDecimal(0, PRECISION);
				}
			}
		}

		for (int a = pivotPositions.size() - 1; a >= 0; a--) {
			for (int b = pivotPositions.get(a) - 1; b >= 0; b--) { // Set left of column = 0
				input[a][b] = new BigDecimal(0, PRECISION);
			}

			// Divide row by pivot position
			divideRow(input, a, input[a][pivotPositions.get(a)]);
			// multiplyAddRow to each row above
			for (int b = a - 1; b >= 0; b--) {
				multiplyAddRow(input, a, b, input[b][pivotPositions.get(a)].divide(input[a][pivotPositions.get(a)], PRECISION).multiply(new BigDecimal(-1.0, PRECISION), PRECISION));
			}
		}

		return input;
	}

	private static void swapRow(BigDecimal[][] input, int rowA, int rowB) {
		BigDecimal[] temp = input[rowA];
		input[rowA] = input[rowB];
		input[rowB]= temp;
	}

	private static void multiplyAddRow(BigDecimal[][] input, int giver, int receiver, BigDecimal multible) {
		BigDecimal[] temp = new BigDecimal[input[receiver].length];
		for (int a = 0; a < temp.length; a++) {
			temp[a] = input[receiver][a].add(input[giver][a].multiply(multible, PRECISION), PRECISION);
		}
		input[receiver] = temp;
	}

	private static void divideRow(BigDecimal[][] input, int row, BigDecimal divder) {
		for (int a = 0; a < input[row].length; a++) {
			input[row][a] = input[row][a].divide(divder, PRECISION);
		}
	}

	private static BigDecimal[][] transpose(BigDecimal[][] input) {
		BigDecimal[][] output = new BigDecimal[input[0].length][input.length];
		for (int a = 0; a < input.length; a++) {
			for (int b = 0; b < input[a].length; b++) {
				output[b][a] = input[a][b];
			}
		}
		return output;
	}

	private static BigDecimal[][] concat(BigDecimal[][] A, BigDecimal[][] B) {
		if (A.length != B.length) {
			throw new IllegalArgumentException("Matricies do not have some amount of rows!");
		}
		BigDecimal[][] output = new BigDecimal[A.length][A[0].length + B[0].length];
		for (int a = 0; a < output.length; a++) {
			for (int b = 0; b < A[0].length; b++) {
				output[a][b] = A[a][b];
			}
		}
		for (int a = 0; a < output.length; a++) {
			for (int b = A[0].length; b < output[0].length; b++) {
				output[a][b] = B[a][b - A[0].length];
			}
		}
		return output;
	}

	private static BigDecimal[][] matrixMultiply(BigDecimal[][] matrixOne, BigDecimal[][] matrixTwo) {
		BigDecimal[][] output = new BigDecimal[matrixOne.length][matrixTwo[0].length];
		for (int a = 0; a < output.length; a++) {
			for (int b = 0; b < output[a].length; b++) {
				output[a][b] = new BigDecimal(0, PRECISION);
			}
		}
		for (int a = 0; a < matrixTwo[0].length; a++) {
			for (int b = 0; b < matrixOne.length; b++) {
				for (int c = 0; c < matrixTwo.length; c++) {
					output[b][a] = output[b][a].add(matrixOne[b][c].multiply(matrixTwo[c][a], PRECISION), PRECISION);
				}
			}
		}
		return output;
	}
}
