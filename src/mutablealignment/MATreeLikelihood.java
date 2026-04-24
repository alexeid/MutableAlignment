package mutablealignment;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

@Description("Tree likelihood that can efficiently recalculate changes in a mutable alignment")
public class MATreeLikelihood extends TreeLikelihood {
	
	private MutableAlignment alignment;
	private boolean alignmentNeedsUpdate;
	private int[] cachedStates;
	private double[] cachedPartials;

	// Tracks tip nodes whose states were transiently overwritten by
	// getLogProbs*Sequence during a proposal, so restore()/accept() can re-sync
	// them from the (post-store/restore) alignment.
	private final Set<Integer> tempTipNodes = new HashSet<>();

	// Tracks internal nodes whose partials buffer has already been flipped during
	// the current proposal. Each node must be flipped at most once per proposal:
	// flipping twice rotates back to the stored slot and overwrites it.
	private final Set<Integer> flippedNodes = new HashSet<>();

	@Override
	public void initAndValidate() {
		if (!(dataInput.get() instanceof MutableAlignment)) {
			throw new IllegalArgumentException("Expected MutableAlignment as data, not " + dataInput.get().getClass().getName());
		}
		alignment = (MutableAlignment) dataInput.get();

		boolean useJava = System.getProperty("java.only") == null ? false : Boolean.valueOf(System.getProperty("java.only"));
		System.setProperty("java.only", "true");
		super.initAndValidate();
		System.setProperty("java.only", useJava + "");

		alignmentNeedsUpdate = false;

		int patternCount = alignment.getPatternCount();
		int stateCount = alignment.getDataType().getStateCount();
		cachedStates = new int[patternCount];
		cachedPartials = new double[patternCount * stateCount];
	}
		
	@Override
	public double calculateLogP() {
		if (alignmentNeedsUpdate) {
			updateAlignment();
			alignmentNeedsUpdate = false;
		}
		// Undo probe-time partials toggles before super.calculateLogP() runs.
		// Reason: traverse() calls setNodePartialsForUpdate() for every dirty
		// ancestor, which toggles the partials-index again. Without this undo,
		// ancestors toggled by getLogProbs*Sequence get toggled twice in one
		// proposal, rotating current back to the stored slot and clobbering it.
		for (Integer n : flippedNodes) {
			likelihoodCore.setNodePartialsForUpdate(n);
		}
		flippedNodes.clear();
		logP = super.calculateLogP();
		
		
//		int [] states = new int[5];
//		for (int i = 0; i < 5; i++) {
//			likelihoodCore.getNodeStates(i, states);
//			System.out.println(i + ": " + Arrays.toString(states));
//		}
		
		return logP;
	}
	
	
	private List<Integer> dirtySequences = new ArrayList<>();
	
	private void updateAlignment() {
        dirtySequences.clear();
        for (Integer i : alignment.getDirtySequenceIndices()) {
        	dirtySequences.add(i);
        }
        updateTipData();
	}

	private void updateTipData() {
        TreeInterface tree = treeInput.get();
    	int patternCount = alignment.getPatternCount();
        int stateCount = alignment.getDataType().getStateCount();
    	for (int nodeNr: dirtySequences) {
    		Node node = tree.getNode(nodeNr);
            int taxonIndex = alignment.getTaxonIndex(node.getID());

            if (m_useAmbiguities.get()) {
	            int k = 0;
	            for (int patternIndex_ = 0; patternIndex_ < patternCount; patternIndex_++) {
	                double[] tipLikelihoods = alignment.getTipLikelihoods(taxonIndex, patternIndex_);
	                if (tipLikelihoods != null) {
	                	for (int state = 0; state < stateCount; state++) {
	                		cachedPartials[k++] = tipLikelihoods[state];
	                	}
	                } else {
	                	int statex = alignment.getPattern(taxonIndex, patternIndex_);
		                boolean[] stateSet = alignment.getStateSet(statex);
		                for (int state = 0; state < stateCount; state++) {
		                	 cachedPartials[k++] = (stateSet[state] ? 1.0 : 0.0);
		                }
	                }
	            }
	            likelihoodCore.setNodePartials(nodeNr, cachedPartials);
	            
            } else {
                DataType dataType = alignment.getDataType();
                for (int i = 0; i < patternCount; i++) {
                    int code = alignment.getPattern(taxonIndex, i);
                    int[] statesForCode = dataType.getStatesForCode(code);
                    if (statesForCode.length==1)
                        cachedStates[i] = statesForCode[0];
                    else
                        cachedStates[i] = code; // Causes ambiguous states to be ignored.
                }
                likelihoodCore.setNodeStates(nodeNr, cachedStates);
                node.makeDirty(Tree.IS_DIRTY);
            }
        }
    }

	/*
	 * returns pattern log likelihoods after setting sequence for node with given nodeNr 
	 * to states encoded in sites
	 */
	public double [] getLogProbsForStateSequence(int nodeNr, int [] sites) {
		// update data for node
		int patternCount = sites.length;
        if (m_useAmbiguities.get()) {
        	throw new IllegalArgumentException("should not use getLogProbsForSequence but getLogProbsForPartialsSequence instead");
        }
        DataType dataType = alignment.getDataType();
        for (int i = 0; i < patternCount; i++) {
            int code = sites[i];
            int[] statesForCode = dataType.getStatesForCode(code);
            if (statesForCode.length==1)
                cachedStates[i] = statesForCode[0];
            else
                cachedStates[i] = code; // Causes ambiguous states to be ignored.
        }
        likelihoodCore.setNodeStates(nodeNr, cachedStates);
        tempTipNodes.add(nodeNr);

        return calcPatternLogLikelihoods(nodeNr);
	}

	/*
	 * returns pattern log likelihoods after setting sequence for node with given nodeNr
	 * to states encoded in sites
	 */
	public double [] getLogProbsForPartialsSequence(int nodeNr, double [] tipLikelihoods) {
        likelihoodCore.setNodePartials(nodeNr, tipLikelihoods);
        tempTipNodes.add(nodeNr);

        return calcPatternLogLikelihoods(nodeNr);
	}

	
	// propagate changes from a leaf node set by getLogProbsForStateSequence or 
	// getLogProbsForPartialsSequence to the root and return updated pattern log 
	// likelihoods
	private double [] calcPatternLogLikelihoods(int nodeNr) {
        // calculate partials up to the root
        Node node = treeInput.get().getNode(nodeNr);
        do {
        	node = node.getParent();
        	final int parentNr = node.getNr();
        	// Flip to scratch buffer so we don't overwrite partials captured by store().
        	// Flip at most once per node per proposal (toggling twice would clobber the stored slot).
        	if (flippedNodes.add(parentNr)) {
        		likelihoodCore.setNodePartialsForUpdate(parentNr);
        	}
            likelihoodCore.calculatePartials(node.getLeft().getNr(), node.getRight().getNr(), parentNr);
        } while (!node.isRoot());
        
        // do fiddly bits at the root
        final double[] proportions = m_siteModel.getCategoryProportions(node);
        likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

        if (getConstantPattern() != null) { // && !SiteModel.g_bUseOriginal) {
            proportionInvariant = m_siteModel.getProportionInvariant();
            // some portion of sites is invariant, so adjust root partials for this
            for (final int i : getConstantPattern()) {
                m_fRootPartials[i] += proportionInvariant;
            }
        }

        // combine with root frequencies
        double[] rootFrequencies = substitutionModel.getFrequencies();
        if (rootFrequenciesInput.get() != null) {
            rootFrequencies = rootFrequenciesInput.get().getFreqs();
        }
        likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);

		return getPatternLogLikelihoods();
	}
	
	
	@Override
	public void store() {
    	dirtySequences.clear();
    	tempTipNodes.clear();
    	flippedNodes.clear();

    	super.store();
	}

	@Override
	public void restore() {
		super.restore();

    	// Resync every tip we temporarily mutated (via getLogProbs*Sequence) from
    	// the alignment, which has itself just been rolled back. This is a superset
    	// of alignment.getDirtySequenceIndices() because getLogProbs* may probe tips
    	// that weren't alignment-dirty (e.g. the partials-fixup leaf).
    	dirtySequences.clear();
    	dirtySequences.addAll(tempTipNodes);
    	for (Integer i : alignment.getDirtySequenceIndices()) {
    		dirtySequences.add(i);
    	}
    	updateTipData();
    	dirtySequences.clear();
    	tempTipNodes.clear();
    	flippedNodes.clear();
	}

	@Override
	protected void accept() {
    	dirtySequences.clear();
    	tempTipNodes.clear();
    	flippedNodes.clear();
		alignment.accept();
		super.accept();
	}

	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty =  super.requiresRecalculation();
		if (alignment.somethingIsDirty()) {
			alignmentNeedsUpdate = true;
            hasDirt = Tree.IS_DIRTY;
			isDirty = true;
		}
		return isDirty;
	}
	
}
