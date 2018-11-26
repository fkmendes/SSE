package SSE;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.google.common.base.Joiner;
import beast.core.Input;
import beast.core.Input.Validate;

/*
 * Note: getCell is where dirty check is done
 */
public class HiddenInstantaneousRateMatrix extends InstantaneousRateMatrix {

	final public Input<Integer> NHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<Boolean> DisallowDoubleTransitionsInput = new Input<>("disallowDoubleTransitions", "Whether or not to set double transition parameters to zero.", Validate.REQUIRED);
	final public Input<Integer> SymmetrifyAcrossDiagonalInput = new Input<>("symmetrifyAcrossDiagonal", "Whether or not to set transition from observed to hidden the same as transitions from hidden to observed.", Validate.REQUIRED);
	final public Input<HiddenObservedStateMapper> HiddenObservedStateMapperInput = new Input<>("hiddenObsStateMapper", "Maps hidden states onto observed states and vice-versa.");
	
	private int numberOfObsStates; // number of observed states (never changes)
	private int numberOfHiddenStates;
	private boolean disallowDoubleTransitions = true;
	private int symmetrifyAcrossDiagonalStateIdx;
	private HiddenObservedStateMapper hiddenObsStateMapper;
	private Map<Integer, List<Integer>> doubleTransitionIRMCellsMap;
	private Double[] matrixContent; // these transition rates will be used to populate 'q' 2D-array
	private Map<Integer, int[]> realParameterToQCell;
	
	@Override
	public void initAndValidate() {
		numberOfObsStates = NstatesInput.get();
		numberOfHiddenStates = NHiddenStatesInput.get(); // when rjMCMC implemented, this can change at different steps
		numberOfStates = numberOfHiddenStates + numberOfObsStates; // this will also change as a result
				
		if (numberOfHiddenStates > 0 && HiddenObservedStateMapperInput.get() == null) {
			throw new IllegalArgumentException("Number of hidden states > 0, but no mapping between observed and hidden states was found.");
		}
		
		disallowDoubleTransitions = DisallowDoubleTransitionsInput.get();
		if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
			hiddenObsStateMapper = HiddenObservedStateMapperInput.get();
		}

		symmetrifyAcrossDiagonalStateIdx = SymmetrifyAcrossDiagonalInput.get(); // if -1, no symmetrifying done
		System.out.println("My symmetrifyAcrossDiagonalStateIdx is " + symmetrifyAcrossDiagonalStateIdx);
		realParameterToQCell = new HashMap<>();
		
		Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
      
		this.populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfObsStates, numberOfHiddenStates, someMatrixContent);
	}
	
	public void populateIRM(boolean ignoreDiagonal, boolean disallowDoubleTransitions, int symmetrifyAcrossDiagonal, int nObsStates, int nHiddenStates, Double[] someMatrixContent) {
		matrixContent = someMatrixContent; 
		numberOfHiddenStates = nHiddenStates;
		numberOfStates = nObsStates + nHiddenStates; // this will also change as a result
		q = new Double[numberOfStates][numberOfStates]; // and so we need to vary the size of q	
		int numberOfInputElements = FlatQmatrixInput.get().getDimension();
		int diagEntry = 0;
		
		// Making sure the right number of transition parameters was provided
		if (ignoreDiagonal) {
			if (!disallowDoubleTransitions) {
				String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements, but not double transitions), but was given ";
				if (((numberOfStates) * (numberOfStates) - (numberOfStates)) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values.");
				} else if (((numberOfStates) * (numberOfStates) - (numberOfStates)) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
			
			else if ((numberOfHiddenStates > 0) && (numberOfHiddenStates == numberOfObsStates) && disallowDoubleTransitions){
				String errorMsg = "Tried to fill transition matrix Q with Nobs=Nhidden (ignoring diagonal elements and also double transitions), but was given ";
				if (((numberOfStates) * (numberOfStates) - (numberOfStates) - (2 * (numberOfObsStates * numberOfHiddenStates - numberOfHiddenStates))) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values."); // subtract diagonal and then double transitions (twice b/c it's top-right and bottom-left)
				} else if (((numberOfStates) * (numberOfStates) - (numberOfStates) - (2 * (numberOfObsStates * numberOfHiddenStates - numberOfHiddenStates))) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
		}
				
		// grabbing cells with double transitions (populating doubleTransitionIRMCellsMap), which will be ignored when populating IRM and remain 0's
		if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
			doubleTransitionIRMCellsMap = new HashMap<>();
			findDoubleTransitionIRMCells(numberOfHiddenStates, numberOfStates);
			// printDoubleTransitionIRMCellsMap(); // for debugging
		}
		
		int realParameterIdx = 0;
		List<Integer> hiddenIdxList = new ArrayList<Integer>(); // will store list of hidden idx associated to observed idx (and be rewritten in a loop)	
		for (int i=0; i<numberOfStates; ++i) {
			for (int j=0; j<numberOfStates; ++j) {
				if (ignoreDiagonal && j == diagEntry) { 
					setCell(i, j, 0.0);
					continue;
				} // ignore diagonal element (make it zero)
				
				if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
					hiddenIdxList = doubleTransitionIRMCellsMap.get(i);
					
					if (hiddenIdxList != null) {
						if (hiddenIdxList.contains(j)) {
							setCell(i, j, 0.0);
							continue;
						}
					}
				} // ignore double transitions (make it zero)
				
				setCell(i, j, someMatrixContent[realParameterIdx]); // after setting zeros (above), we finally set the matrix cells
				
				/* 
				 * recording which cells of the Q matrix we put some non-zero value in (mapping cell to index in RealParameter list
				 * so we can use it later to resize matrixContent in rjMCMC)
				 */
				if (realParameterToQCell.get(realParameterIdx) == null) {
					int[] nonZeroIJ = {i, j};
					realParameterToQCell.put(realParameterIdx, nonZeroIJ);
					// System.out.println("HIRM: recording RealParameter transition index " + realParameterIdx + " and its cell " + Arrays.toString(nonZeroIJ));
				}
				
				realParameterIdx += 1;
			}
				 
			diagEntry += 1; // index of the diagonal element we should ignore
		}
		
		// make q_hidden2obs = q_obs2hidden for all states
		if (symmetrifyAcrossDiagonal != -1) {
			symmetrifyAcrossDiagonal(symmetrifyAcrossDiagonal);
		}
		
		irmDirty = false; // after re-population of IRM, things are clean
	}
	
	/* Later I should make sure this is made just once for each k when rjMCMC is implemented, or even when any transition rate is operated on */
	public void findDoubleTransitionIRMCells(int numberOfHiddenStates, int totalNumberOfStates) {
		
		for (int i=0; i<totalNumberOfStates; ++i) {
			// i-th row

			// recording double transition cells from top-right of IRM
			if (i < numberOfObsStates) {
				for (int j=numberOfObsStates; j<totalNumberOfStates; ++j) {
					// j-th col

					int h = hiddenObsStateMapper.getHiddenFromObs(i); // hidden state idx
					
					// this is checking the diagonal of the hidden states sub-matrix (off-diagonal are 0!)
					// if double transition; or if Nhidden < Nobs, and this obs never transitions into a hidden state (h == -1)
					if ((j - numberOfObsStates) != h || h == -1) {
					
						// j is hidden state idx
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
					}
					
				}
			}
			
			// recording double transition cells from bottom-left of IRM
			if (i >= numberOfObsStates) {
				int h = i-numberOfObsStates; // hidden state idx
				Collection<Integer> obsIdxColl = hiddenObsStateMapper.getObsCollFromHidden(h); // collection of observed state idxs
				
				// left cols (j) of the bottom half
				for (int j=0; j < numberOfObsStates; ++j) {
					/* the observed states that never transition to/from hidden states are taken care of automatically below,
					 * as their indexes will not be obsIdxColl
					 */
					if (!obsIdxColl.contains(j)) {		
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
					}

					// System.out.println("Bottom-left i=" + i + " j=" +j);
				}
			}
		}
					
	}
	
	/* necessary for the intermediate step in rjMCMC */
	public void symmetrifyAcrossDiagonal(int obsStateToSymmIdx) {

		System.out.println("Symmetrifying observed/hidden state " + obsStateToSymmIdx);
		
		for (int col = numberOfObsStates; col < numberOfStates; ++col) {
			if (q[obsStateToSymmIdx][col] > 0.0) {
				setCell(col, obsStateToSymmIdx, q[obsStateToSymmIdx][col]);
			}
		}
	}
	
	public void resetQ(int aNumberOfStates, int aNumberOfHiddenStatesInMask) {
		numberOfStates = aNumberOfStates + aNumberOfHiddenStatesInMask;
		q = new Double[numberOfStates][numberOfStates];
	}

	/* 
	 * called by masquerade ball during rjMCMC
	 */
	public void setHiddenObsStateAssignment(int[] aHiddenToObsStateAssignment) {
		hiddenObsStateMapper.setHiddenStateStrings(aHiddenToObsStateAssignment);
	}
	
	// getters
	@Override
	public int getNumObsStates() {
		return numberOfObsStates;
	}
	
	public int getNumHiddenStates() {
		return numberOfHiddenStates;
	}
	
	public int getNumStates () {
		return numberOfStates; // total number of states
	}
	
	@Override
	public double getCell(int from, int to, double rate) {
		synchronized (this) {
			if (irmDirty) {
				Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
				populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfObsStates, numberOfHiddenStates, someMatrixContent); // only re-populate IRM if some transition rate was operated on
			}			
		}
		return q[from][to] * rate;
	}

	public Double[] getMatrixContent() {
		return matrixContent;
	}
	
	public Map<Integer, int[]> getRealParameterToQCellMap() {
		return realParameterToQCell;
	}
	
	/*
	 * only used by masquerade ball (not by ODE),
	 * so does not check for dirtyness nor populateIRM
	 * (getCell is the method that does that because that's what the ODE uses)
	 */
	public Double[][] getQ() {
		return q;
	}
	
	public void printDoubleTransitionIRMCellsMap() {
		System.out.println("IRM cells that contain double transitions and that should get 0s.");
		Iterator<Entry<Integer, List<Integer>>> iter = doubleTransitionIRMCellsMap.entrySet().iterator();
		while (iter.hasNext()) {
			Entry<Integer, List<Integer>> entry = iter.next();
			String k = Integer.toString(entry.getKey());
			String v = Joiner.on(",").join(entry.getValue());
			System.out.println("i=" + k + " j(s)=" + v);
		}
	}
	
	public void printHiddenObservedStateMapper() {
		hiddenObsStateMapper.printHidden2ObsMap();
		hiddenObsStateMapper.printObs2HiddenMap();
	}

	@Override
	protected void restore() {
		Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
		populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfObsStates, numberOfHiddenStates, someMatrixContent);
		super.restore();
	}
}