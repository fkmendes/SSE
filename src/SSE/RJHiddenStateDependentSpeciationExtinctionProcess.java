package SSE;

import java.util.Arrays;
import java.util.concurrent.Executors;

import beast.app.BeastMCMC;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Tree;

public class RJHiddenStateDependentSpeciationExtinctionProcess extends HiddenStateDependentSpeciationExtinctionProcess {

	final public Input<MasqueradeBall> masqueradeBallInput = new Input<>("masqueradeBall", "MaskeradeBall object that allows resetting lambdas, mus and Q's depending on the model mask.", Validate.REQUIRED);
	
	private MasqueradeBall masqueradeBall;
	
	public RJHiddenStateDependentSpeciationExtinctionProcess() {
		hirmInput.setRule(Input.Validate.FORBIDDEN);
		lambdaMuAssignerInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		masqueradeBall = masqueradeBallInput.get();
		tree = treeInput.get();
		hiddenTraitStash = hiddenTraitStashInput.get();
		incorporateCladogenesis = cladoFlagInput.get();
		
		pi = masqueradeBall.getPis();
		mu = masqueradeBall.getMus();
		q = masqueradeBall.getHIRM();
		
		// likelihood-related
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		scalingConstants = new double[tree.getNodeCount()]; // equivalent to diversitree's lq (but not in log-scale), these are used as denominators during the likelihood computation
		Arrays.fill(scalingConstants, 1.0);		
		finalLk = 0.0;
		finalLogLk = 0.0;
		
		// cache-related 
		storedNodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2];  
		storedScalingConstants = new double[tree.getNodeCount()]; 
		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		hasDirt = Tree.IS_FILTHY;

		useThreads = useThreadsInput.get() && (BeastMCMC.m_nThreads > 1);
		nrOfThreads = useThreads ? BeastMCMC.m_nThreads : 1;
		if (useThreads && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		}
		if (useThreads) {
			exec = Executors.newFixedThreadPool(nrOfThreads);
		}
	}
	
	@Override
	public double calculateLogP() {
		pi = masqueradeBall.getPis();
		mu = masqueradeBall.getMus();
		q = masqueradeBall.getHIRM();

		// when a character state parameter changes, tree is filthy, we can use threads
		// otherwise, caching allows us to only recompute part of the tree likelihood, and thread overhead not worth it
		if (hasDirt == Tree.IS_FILTHY && useThreads) {
			computeNodeLkUsingThreads();
		} else {
			computeNodeLk(tree.getRoot(), true);
		}
		
		logP = finalLogLk;
		return logP;
	}
	
}