package SSE;

public class HiddenSSEODE extends SSEODE {

	private HiddenInstantaneousRateMatrix hq; // ctor arg
	
	/*
	 * Constructor
	 * Speciation rates and event map are set independently so more or less general models can use this class
	 */
	public HiddenSSEODE(Double[] mu, HiddenInstantaneousRateMatrix hq, double rate, boolean incorporateCladogenesis) {
		super(mu, rate, incorporateCladogenesis);
		this.hq = hq;
		numStates = hq.getNumStates();
	}
	
	@Override
	protected double getQCell(int from, int to, double rate) {
		return hq.getCell(from, to, rate);
	}
}