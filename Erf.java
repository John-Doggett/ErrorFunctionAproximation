import java.math.BigDecimal;
import java.math.MathContext;

public class Erf implements Function{

	@Override
	public BigDecimal evaulate(BigDecimal X, int precision) {
		MathContext PRECISION = new MathContext(precision);
		BigDecimal sum = new BigDecimal(0, PRECISION);
		BigDecimal n = new BigDecimal(5000,PRECISION);
		BigDecimal increment = (X).divide(n, PRECISION);
		BigDecimal four = new BigDecimal(4, PRECISION);
		BigDecimal two = new BigDecimal(2, PRECISION);

		sum = sum.add(new BigDecimal(Math.exp(0), PRECISION));
		for(double a = 1.0; a < n.doubleValue(); a = a + 1.0) {
			if(a%2.0==1.0) {
				sum = sum.add(four.multiply(new BigDecimal(Math.exp(-Math.pow(a*increment.doubleValue(), 2)), PRECISION), PRECISION), PRECISION);
			}
			else {
				sum = sum.add(two.multiply(new BigDecimal(Math.exp(-Math.pow(a*increment.doubleValue(), 2)), PRECISION), PRECISION), PRECISION);
			}
		}
		sum = sum.add(new BigDecimal(Math.exp(-Math.pow(X.doubleValue(), 2)), PRECISION), PRECISION);

		sum = sum.multiply((X).divide(n, PRECISION),PRECISION);
		sum = sum.divide(new BigDecimal(3.0,PRECISION),PRECISION);
		sum = sum.multiply( two.divide(new BigDecimal(Math.PI).sqrt(PRECISION), PRECISION), PRECISION);
		return sum;
	}

}
