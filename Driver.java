import java.math.BigDecimal;

public class Driver{

	public static void main(String[] args) {
		Function fun = new Erf();
		Long time = System.nanoTime();
		System.out.println(fun.evaulate(new BigDecimal(1), 50).doubleValue());
		System.out.println("Took: " + (System.nanoTime()-time));
		double[] polynomial = Aproximate.aproximate(fun,0,4,100, 40, 50);
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
