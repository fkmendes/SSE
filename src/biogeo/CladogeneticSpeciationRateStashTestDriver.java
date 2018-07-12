package biogeo;

import java.util.ArrayList;
import java.util.List;

import beast.core.parameter.RealParameter;
import biogeo.CladoTriplet.speciationType;

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
		cladoTripletList.add(sTriplet);
		cladoTripletList.add(ssTriplet);
		
		Double[] sSpeciationRate = {0.1};
		Double[] ssSpeciationRate = {0.2};
		Double[] vSpeciationRate = {0.0};
		RealParameter sympatricSpeciationRate = new RealParameter(sSpeciationRate);
		RealParameter subSympatricSpeciationRate = new RealParameter(ssSpeciationRate);
		RealParameter vicariantSpeciationRate = new RealParameter(vSpeciationRate);

		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("CladoTriplets", cladoTripletList,
				"SympatricRate", sympatricSpeciationRate,
				"SubsympatricRate", subSympatricSpeciationRate,
				"VicariantRate", vicariantSpeciationRate);
		csrt.printEventMap();
		
//		int[][] cladogenetic_events = {{2, 2, 1},{1, 1, 1}};
//		double[] lambda = {0.1, 0.2};
//		CladogeneticSpeciationRateStash clado_stash = new CladogeneticSpeciationRateStash(cladogenetic_events, lambda);
//		clado_stash.printEventMap();		
	}
}

