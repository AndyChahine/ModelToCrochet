
public class Params {
	int totalEpochs;
	double sumOfErrors;
	double meanOfErrors;
	double c;
	double rho;
	
	public void update(int totalEpochs, double sumOfErrors, double meanOfErrors, double c, double rho) {
		this.totalEpochs = totalEpochs;
		this.sumOfErrors = sumOfErrors;
		this.meanOfErrors = meanOfErrors;
		this.c = c;
		this.rho = rho;
	}
	
	@Override
	public String toString() {
		return "epochs: " + totalEpochs + " sum of errors: " + sumOfErrors + " mean error: " + meanOfErrors + " c: " + c + " rho: " + rho;
	}
}
