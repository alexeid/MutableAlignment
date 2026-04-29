package mutablealignment;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import beagle.Beagle;
import beast.base.core.Description;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;

@Description("Tree likelihood that can efficiently recalculate changes in a mutable alignment")
public class BeagleMATreeLikelihood extends BeagleTreeLikelihood {
	
	private MutableAlignment alignment;
	private boolean alignmentNeedsUpdate;
	private int[] cachedOperations;
	private int[] cachedStates;
	private double[] cachedPartials;

	// Mapping from alignment column index to tree leaf node number, and
	// inverse. Computed once in initAndValidate by taxon name; alignment
	// column order and tree leaf order are independent.
	private int[] alignmentIdxToTreeNodeNr;
	private int[] treeNodeNrToAlignmentIdx;

	// Tracks tip nodes whose states were transiently overwritten by
	// getLogProbs*Sequence during a proposal, so restore()/accept() can re-sync
	// them from the (post-store/restore) alignment. BEAGLE's tip states are
	// not double-buffered, so this resync mechanism is unavoidable.
	// Stored as alignment column indices, matching dirtySequences.
	private final Set<Integer> tempTipNodes = new HashSet<>();

	@Override
	public void initAndValidate() {
		if (!(dataInput.get() instanceof MutableAlignment)) {
			throw new IllegalArgumentException("Expected MutableAlignment as data, not " + dataInput.get().getClass().getName());
		}
		alignment = (MutableAlignment) dataInput.get();

		super.initAndValidate();

		alignmentNeedsUpdate = false;

		int patternCount = alignment.getPatternCount();
		int stateCount = alignment.getDataType().getStateCount();
		cachedStates = new int[patternCount];
		cachedPartials = new double[patternCount * stateCount];
		cachedOperations = new int[treeInput.get().getNodeCount() * Beagle.OPERATION_TUPLE_SIZE];

		buildTaxonIndexMaps();
	}

	private void buildTaxonIndexMaps() {
		int taxonCount = alignment.getTaxonCount();
		alignmentIdxToTreeNodeNr = new int[taxonCount];
		treeNodeNrToAlignmentIdx = new int[taxonCount];
		TreeInterface tree = treeInput.get();
		for (Node leaf : tree.getExternalNodes()) {
			int nodeNr = leaf.getNr();
			int alignIdx = alignment.getTaxonIndex(leaf.getID());
			if (alignIdx < 0) {
				throw new RuntimeException("Tree leaf " + leaf.getID() + " not found in alignment");
			}
			alignmentIdxToTreeNodeNr[alignIdx] = nodeNr;
			treeNodeNrToAlignmentIdx[nodeNr] = alignIdx;
		}
	}
		
	@Override
	public double calculateLogP() {
		if (alignmentNeedsUpdate) {
			updateAlignment();
			alignmentNeedsUpdate = false;
		}
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
    	for (int taxonIndex: dirtySequences) {
    		// dirtySequences holds alignment column indices.
    		int nodeNr = alignmentIdxToTreeNodeNr[taxonIndex];
    		Node node = tree.getNode(nodeNr);

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
	            beagle.setPartials(nodeNr, cachedPartials);
	            
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
                beagle.setTipStates(nodeNr, cachedStates);
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
        beagle.setTipStates(nodeNr, cachedStates);
        tempTipNodes.add(treeNodeNrToAlignmentIdx[nodeNr]);

        return calcPatternLogLikelihoods(nodeNr);
	}

	/*
	 * returns pattern log likelihoods after setting sequence for node with given nodeNr
	 * to states encoded in sites
	 */
	public double [] getLogProbsForPartialsSequence(int nodeNr, double [] tipLikelihoods) {
        beagle.setPartials(nodeNr, tipLikelihoods);
        tempTipNodes.add(treeNodeNrToAlignmentIdx[nodeNr]);

        return calcPatternLogLikelihoods(nodeNr);
	}

	
	// propagate changes from a leaf node set by getLogProbsForStateSequence or
	// getLogProbsForPartialsSequence to the root and return updated pattern log
	// likelihoods. Hermetic: flips each touched ancestor's partials offset to
	// the scratch slot before writing, computes pattern log likelihoods at the
	// root, then flips each ancestor back. After this method returns the
	// partialBufferHelper offsets are exactly what they were on entry.
	private double [] calcPatternLogLikelihoods(int nodeNr) {
		return calcPatternLogLikelihoods(nodeNr, new HashSet<>());
	}

	// Inner overload threads the per-call `flipped` set through the rescale
	// recursion so we still flip each ancestor at most once across attempts
	// (and flip them all back on the final success path).
	private double [] calcPatternLogLikelihoods(int nodeNr, Set<Integer> flipped) {

        Node node = treeInput.get().getNode(nodeNr);
        int operationCount = 0;
        final int[] operations = cachedOperations;
        do {
        	node = node.getParent();
        	nodeNr = node.getNr();

            // Traverse down the two child nodes
            Node child1 = node.getLeft();

            Node child2 = node.getRight();

            int x = operationCount * Beagle.OPERATION_TUPLE_SIZE;

            // Flip to scratch slot so we don't overwrite partials captured by store().
            // Flip at most once per node per proposal (toggling twice would clobber the stored slot).
            if (flipped.add(nodeNr)) {
                partialBufferHelper.flipOffset(nodeNr);
            }
            operations[x] = partialBufferHelper.getOffsetIndex(nodeNr);

            if (useScaleFactors) {
                // get the index of this scaling buffer
                int n = nodeNr - tipCount;

                if (recomputeScaleFactors) {
                    // flip the indicator: can take either n or (internalNodeCount + 1) - n
                    // scaleBufferHelper.flipOffset(n);

                    // store the index
                    scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                    operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                    operations[x + 2] = Beagle.NONE;

                } else {
                    operations[x + 1] = Beagle.NONE;
                    operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                }

            } else {

                if (useAutoScaling) {
                    scaleBufferIndices[nodeNr - tipCount] = partialBufferHelper.getOffsetIndex(nodeNr);
                }
                operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                operations[x + 2] = Beagle.NONE;
            }

            operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNr()); // source node 1
            operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNr()); // source matrix 1
            operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNr()); // source node 2
            operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNr()); // source matrix 2

            operationCount++;

		} while (!node.isRoot());
        
        
        
        
        double logL;
        boolean done;
        boolean firstRescaleAttempt = true;

        do {

            beagle.updatePartials(operations, operationCount, Beagle.NONE);

            Node root = treeInput.get().getRoot();
            
            int rootIndex = partialBufferHelper.getOffsetIndex(root.getNr());

            double[] categoryWeights = m_siteModel.getCategoryProportions(null);
            if (getConstantPattern() != null) {
	            double [] tmp = new double [categoryWeights.length - 1];
	            for (int k = 0; k < invariantCategory; k++) {
	            	tmp[k] = categoryWeights[k];
	            }
	            for (int k = invariantCategory + 1; k < categoryWeights.length; k++) {
	            	tmp[k-1] = categoryWeights[k];
	            }
	            categoryWeights = tmp;
            }
            double[] frequencies = rootFrequenciesInput.get() == null ?
                    				substitutionModel.getFrequencies() :
                    				rootFrequenciesInput.get().getFreqs();


            int cumulateScaleBufferIndex = Beagle.NONE;
            if (useScaleFactors) {

                if (recomputeScaleFactors) {
                    scaleBufferHelper.flipOffset(internalNodeCount);
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                    beagle.resetScaleFactors(cumulateScaleBufferIndex);
                    beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, cumulateScaleBufferIndex);
                } else {
                    cumulateScaleBufferIndex = scaleBufferHelper.getOffsetIndex(internalNodeCount);
                }
            } else if (useAutoScaling) {
                beagle.accumulateScaleFactors(scaleBufferIndices, internalNodeCount, Beagle.NONE);
            }

            // these could be set only when they change but store/restore would need to be considered
            
            for (int i = 0; i < categoryWeights.length; i++) {
            	if (categoryWeights[i] != currentCategoryWeights[i]) {
                    beagle.setCategoryWeights(0, categoryWeights);
            		i = categoryWeights.length;
            	}
            }
            currentCategoryWeights = categoryWeights;
            for (int i = 0; i < frequencies.length; i++) {
            	if (frequencies[i] != currentFreqs[i]) {
                    beagle.setStateFrequencies(0, frequencies);
            		i = frequencies.length;
            	}
            }
            currentFreqs = frequencies;

            double[] sumLogLikelihoods = new double[1];

            beagle.calculateRootLogLikelihoods(new int[]{rootIndex}, new int[]{0}, new int[]{0},
                    new int[]{cumulateScaleBufferIndex}, 1, sumLogLikelihoods);

            logL = sumLogLikelihoods[0];

            beagle.getSiteLogLikelihoods(patternLogLikelihoods);
            if (ascertainedSitePatterns) {
                // Need to correct for ascertainedSitePatterns
                logL = getAscertainmentCorrectedLogLikelihood(alignment,
                        patternLogLikelihoods, alignment.getWeights(), frequencies);
            } else if (invariantCategory >= 0) {
                int [] patternWeights = alignment.getWeights();
                proportionInvariant = m_siteModel.getProportionInvariant();
                
                
    	        for (int k : getConstantPattern()) {
    	        	int i = k / m_nStateCount;
    	        	int j = k % m_nStateCount;
    	        	patternLogLikelihoods[i] = (Math.log(Math.exp(patternLogLikelihoods[i]) + proportionInvariant * frequencies[j]));
    	        }
        	
	            logL = 0.0;
	            for (int i = 0; i < patternCount; i++) {
	                logL += patternLogLikelihoods[i] * patternWeights[i];
	            }
            }

            if (Double.isNaN(logL) || Double.isInfinite(logL)) {
                everUnderflowed = true;
                logL = Double.NEGATIVE_INFINITY;

                if (firstRescaleAttempt && (rescalingScheme == PartialsRescalingScheme.DYNAMIC || rescalingScheme == PartialsRescalingScheme.DELAYED)) {
                    // we have had a potential under/over flow so attempt a rescaling                	
                	useScaleFactors = true;
                    recomputeScaleFactors = true;

                    for (int i = 0; i < eigenCount; i++) {
                        branchUpdateCount[i] = 0;
                    }

                    operationCount = 0;

                    // traverse again but without flipping partials indices as we
                    // just want to overwrite the last attempt. We will flip the
                    // scale buffer indices though as we are recomputing them.
                    // Pass `flipped` so the recursion's add() check skips
                    // re-flipping the ancestors we already toggled.
                    return calcPatternLogLikelihoods(nodeNr, flipped);

//                    done = false; // Run through do-while loop again
//                    firstRescaleAttempt = false; // Only try to rescale once
                } else {
                    // we have already tried a rescale, not rescaling or always rescaling
                    // so just return the likelihood...
                    done = true;
                }
            } else {
                done = true; // No under-/over-flow, then done
            }

        } while (!done);

        // flip ancestors back so the partials offsets are unchanged on exit
        for (Integer nr : flipped) {
            partialBufferHelper.flipOffset(nr);
        }

        return patternLogLikelihoods.clone();
    }

	
	@Override
	public void store() {
    	dirtySequences.clear();
    	// Do NOT clear tempTipNodes here. store() is called by MCMC between
    	// operator.proposal() and calculateLogP() (default
    	// requiresStateInitialisation=true). Tip-state probes happen during
    	// proposal, and the resync needs to know which tips were touched until
    	// restore()/accept() handles them.

    	super.store();
	}

	@Override
	public void restore() {
		super.restore();

    	// Resync every tip we temporarily mutated from the (now-rolled-back)
    	// alignment. dirtySequences was populated during the previous
    	// calculateLogP from alignment.getDirtySequenceIndices(); we can't
    	// re-query that here because alignment.restore() has already cleared
    	// the edit list. Add probe-touched tips (tempTipNodes) since they may
    	// not have been alignment-dirty (e.g. ExchangeGibbsOperator's partials-
    	// fixup leaf).
    	dirtySequences.addAll(tempTipNodes);
    	updateTipData();
    	dirtySequences.clear();
    	tempTipNodes.clear();
	}

	@Override
	protected void accept() {
    	dirtySequences.clear();
    	tempTipNodes.clear();
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
