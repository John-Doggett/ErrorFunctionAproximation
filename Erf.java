
public class Erf implements Function{

	@Override
	public double evaulate(double X) {
		double sum = 0;
		double n = 3000;
		double increment = (X)/n;
		sum += Math.exp(-Math.pow(0, 2));
		for(double a = 1.0; a < n; a = a + 1.0) {
			if(a%2.0==1.0) {
				sum+=4.0*(Math.exp(-Math.pow(a*increment, 2)));
			}
			else {
				sum+=2.0*(Math.exp(-Math.pow(a*increment, 2)));
			}
		}
		sum += Math.exp(-Math.pow(X, 2));
		sum *= (X)/n;
		sum /= 3.0;
		sum *= 2.0/(Math.sqrt(Math.PI));
		return sum;
	}

}
