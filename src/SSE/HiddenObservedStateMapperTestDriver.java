package SSE;

public class HiddenObservedStateMapperTestDriver {

	public static void main(String[] args) {
		String hiddenStatesString = "0,1,2,2";
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
		System.out.println("hidden from obs=0 is " + stateMapper.getHiddenFromObs(0));
		System.out.println("hidden from obs=3 is " + stateMapper.getHiddenFromObs(3));
		System.out.println("obs from hidden=0 is " + stateMapper.getObsCollFromHidden(0));
		System.out.println("obs from hidden=2 is " + stateMapper.getObsCollFromHidden(2));
	}
}
