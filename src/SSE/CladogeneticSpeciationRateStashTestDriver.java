package SSE;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import SSE.CladoTriplet.speciationType;
import beast.core.parameter.RealParameter;

public class CladogeneticSpeciationRateStashTestDriver {

	public static void main(String[] args) {
		CladoTriplet sTriplet = new CladoTriplet();
		sTriplet.initByName("ParentState", 1,
				"LeftChildState", 1,
				"RightChildState", 1,
				"SpeciationType", speciationType.SYMPATRY);
		
		CladoTriplet ssTriplet = new CladoTriplet();
		ssTriplet.initByName("ParentState", 2,
				"LeftChildState", 2,
				"RightChildState", 1,
				"SpeciationType", speciationType.SUBSYMPATRY);
		
		List<CladoTriplet> cladoTripletList = new ArrayList<CladoTriplet>();
		Collections.addAll(cladoTripletList, sTriplet, ssTriplet);
		
		Double[] sSpeciationRate = {0.1};
		Double[] ssSpeciationRate = {0.2};
		Double[] vSpeciationRate = {0.0};
		Double[] jSpeciationRate = {0.0};
		RealParameter sympatricSpeciationRate = new RealParameter(sSpeciationRate);
		RealParameter subSympatricSpeciationRate = new RealParameter(ssSpeciationRate);
		RealParameter vicariantSpeciationRate = new RealParameter(vSpeciationRate);
		RealParameter jumpSpeciationRate = new RealParameter(jSpeciationRate);

		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("CladoTriplets", cladoTripletList,
				"SympatricRate", sympatricSpeciationRate,
				"SubsympatricRate", subSympatricSpeciationRate,
				"VicariantRate", vicariantSpeciationRate,
				"JumpRate", jumpSpeciationRate);
		csrt.printEventMap();	
	}
}

