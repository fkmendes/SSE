package SSE;

public class HiddenObservedStateMapperTestDriver {

	public static void main(String[] args) {
		String hiddenStatesString = "0,1,2,2";
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		stateMapper.makeMaps();
	}
}
