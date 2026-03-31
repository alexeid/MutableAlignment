package test.mutablealignment;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.Before;
import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.HKY;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.State;
import mutablealignment.MATreeLikelihood;
import mutablealignment.MutableAlignment;
import test.beast.BEASTTestCase;

public class MATreeLikelihoodTest {

	private double logP1, logP2, logP3, logP4, logP5;
	private Tree tree;
	private SiteModel siteModel;
	
	@Before
	public void SetUp() throws Exception {
        logP1 = calcLogP(MutableAlignmentTest.getAlignment1());
        logP2 = calcLogP(MutableAlignmentTest.getAlignment2());
        logP3 = calcLogP(MutableAlignmentTest.getAlignment3());
        logP4 = calcLogP(MutableAlignmentTest.getAlignment4());
        logP5 = calcLogP(MutableAlignmentTest.getAlignment5());        
	}
	
    static public Tree getTree(Alignment data) throws Exception {
        TreeParser tree = new TreeParser();
        tree.initByName("taxa", data,
                "newick", "((((0human:0.024003,(1chimp:0.010772,2bonobo:0.010772):0.013231):0.012035,3gorilla:0.036038):0.033087000000000005,4orangutan:0.069125):0.030456999999999998,5siamang:0.099582);",
                "IsLabelledNewick", true);
        tree.setID("Tree.t:tree");
        return tree;
    }

	private double calcLogP(Alignment data) throws Exception {
        tree = getTree(data);

        Frequencies freqs = new Frequencies();
        freqs.initByName("frequencies", "0.1 0.2 0.3 0.4");

        HKY hky = new HKY();
        hky.initByName("kappa", "20", "frequencies", freqs);

        siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", hky);

        TreeLikelihood likelihood = new TreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);

        double logP = 0;
        logP = likelihood.calculateLogP();
		return logP;
	}
	
	@Test
	public void testMATreeLikelihood() throws Exception {
		MutableAlignment a = MutableAlignmentTest.getAlignment1();
        MATreeLikelihood likelihood = new MATreeLikelihood();
        likelihood.initByName("data", a, "tree", tree, "siteModel", siteModel);
        
        State state = new State();
        state.initByName("stateNode", a);
        state.initialise();
		state.setPosterior(likelihood);
        double logP = state.robustlyCalcPosterior(likelihood);
        assertEquals(logP, logP1, BEASTTestCase.PRECISION);
        
        // make an edit of a single site
        state.store(1);
        state.storeCalculationNodes();
        a.setSiteValue(1, 1, 0);
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP2, logP, BEASTTestCase.PRECISION);
        state.restore();
        state.restoreCalculationNodes();

        // restore/undo edit
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP1, logP, BEASTTestCase.PRECISION);
        state.acceptCalculationNodes();

        // make multiple edits: a single site + taxon value edit 
        state.store(1);
        state.storeCalculationNodes();
        a.setSiteValue(1, 1, 0);
		a.setSiteValuesByTaxon(2, new int[] {3,2,1});
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP3, logP, BEASTTestCase.PRECISION);
        state.restore();
        state.restoreCalculationNodes();

        // restore/undo edit
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP1, logP, BEASTTestCase.PRECISION);
        state.acceptCalculationNodes();

	
        // make multiple edits: a single site + taxon value + site value edit
        state.store(1);
        state.storeCalculationNodes();
        a.setSiteValue(1, 1, 0);
		a.setSiteValuesByTaxon(2, new int[] {3,2,1});
		a.setSiteValuesBySite(0, new int[] {1,1,1,1,1,1});
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP4, logP, BEASTTestCase.PRECISION);
        state.restore();
        state.restoreCalculationNodes();

        // restore/undo edit
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP1, logP, BEASTTestCase.PRECISION);
        state.acceptCalculationNodes();

        // edit whole alignment
        state.store(1);
        state.storeCalculationNodes();
		a.setSiteValues(new int[][]{{2,1,2,1,2,1},{1,2,1,2,1,2},{2,1,2,1,2,1}});
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP5, logP, BEASTTestCase.PRECISION);
        state.restore();
        state.restoreCalculationNodes();

        // restore/undo edit
        state.checkCalculationNodesDirtiness();
        logP = likelihood.calculateLogP();
        assertEquals(logP1, logP, BEASTTestCase.PRECISION);
        state.acceptCalculationNodes();
	}
}
