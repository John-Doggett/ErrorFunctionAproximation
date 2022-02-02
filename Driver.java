
public class Driver{

	public static void main(String[] args) {
		Function fun = new Erf();
		Long time = System.nanoTime();
		System.out.println(fun.evaulate(1));
		System.out.println("Took: " + (System.nanoTime()-time));
		// Lower Bound, Upper Bound, LS Intervals, Polynomial count, Precision (decimal places)
		double[] polynomial = Aproximate.aproximate(fun,0,4,100, 30, 50);
		time = System.nanoTime();
		double x=1.0;
		double sum = 0.0;
		for(int a = 0; a< polynomial.length; a++) {
			sum+=polynomial[a]*Math.pow(x, a);
		}
		time = System.nanoTime()-time;
		System.out.println(sum);
		System.out.println("Took: " + time + "\n");
		
		for(int a = 0; a < polynomial.length; a++) {
			System.out.println(polynomial[a] + " * x^" + a);
		}
	}
}
