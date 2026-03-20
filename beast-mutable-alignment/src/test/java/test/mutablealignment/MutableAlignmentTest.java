package test.mutablealignment;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import mutablealignment.MutableAlignment;

public class MutableAlignmentTest {
	
	@Test
	public void testMutableAlignment() throws Exception {
		MutableAlignment a = getAlignment1();
		a.setSiteValue(1, 1, 0);
		MutableAlignment b = getAlignment2();
		assertAlignmentsAreEqual(a, b);

		a.setSiteValuesByTaxon(2, new int[] {3,2,1});
		MutableAlignment c = getAlignment3();
		assertAlignmentsAreEqual(a, c);

		a.setSiteValuesBySite(0, new int[] {1,1,1,1,1,1});
		MutableAlignment d = getAlignment4();
		assertAlignmentsAreEqual(a, d);

		a.setSiteValues(new int[][]{{2,1,2,1,2,1},{1,2,1,2,1,2},{2,1,2,1,2,1}});
		MutableAlignment e = getAlignment5();
		assertAlignmentsAreEqual(a, e);
	}
	

	static Sequence human = new Sequence("0human",         "AGA");
    static Sequence chimp = new Sequence("1chimp",         "AGA");
    static Sequence bonobo = new Sequence("2bonobo",       "AGA");
    static Sequence gorilla = new Sequence("3gorilla",     "AGA");
    static Sequence orangutan = new Sequence("4orangutan", "AGA");
    static Sequence siamang = new Sequence("5siamang",     "TGA");
	
	
    static public MutableAlignment getAlignment1() throws Exception {
        return getAlignment(human, chimp, bonobo, gorilla, orangutan, siamang);

    }

    static public MutableAlignment getAlignment2() throws Exception {
        Sequence chimp = new Sequence("1chimp",         "AAA"); // one char diff from Alignment1
        return getAlignment(human, chimp, bonobo, gorilla, orangutan, siamang);
    }
	
    static public MutableAlignment getAlignment3() throws Exception {
        Sequence chimp = new Sequence("1chimp",         "AAA"); // one char diff from Alignment1
        Sequence bonobo = new Sequence("2bonobo",       "TGC"); // one sequence diff from Alignment1
        return getAlignment(human, chimp, bonobo, gorilla, orangutan, siamang);
    }

    static public MutableAlignment getAlignment4() throws Exception {
        Sequence human = new Sequence("0human",         "CGA");
        Sequence chimp = new Sequence("1chimp",         "CAA");
        Sequence bonobo = new Sequence("2bonobo",       "CGC");
        Sequence gorilla = new Sequence("3gorilla",     "CGA");
        Sequence orangutan = new Sequence("4orangutan", "CGA");
        Sequence siamang = new Sequence("5siamang",     "CGA");
        return getAlignment(human, chimp, bonobo, gorilla, orangutan, siamang);
    }

    static public MutableAlignment getAlignment5() throws Exception {
        Sequence human = new Sequence("0human",         "GCG");
        Sequence chimp = new Sequence("1chimp",         "CGC");
        Sequence bonobo = new Sequence("2bonobo",       "GCG");
        Sequence gorilla = new Sequence("3gorilla",     "CGC");
        Sequence orangutan = new Sequence("4orangutan", "GCG");
        Sequence siamang = new Sequence("5siamang",     "CGC");
        return getAlignment(human, chimp, bonobo, gorilla, orangutan, siamang);
    }

    private static MutableAlignment getAlignment(Sequence human, Sequence chimp, Sequence bonobo, Sequence gorilla,
			Sequence orangutan, Sequence siamang) {
        MutableAlignment data = new MutableAlignment();
        data.initByName("sequence", human, "sequence", chimp, "sequence", bonobo, "sequence", gorilla, "sequence", orangutan, "sequence", siamang,
                "dataType", "nucleotide"
        );
        return data;
	}

	static void assertAlignmentsAreEqual(Alignment a, Alignment b) {
		assertEquals(a.getPatternCount(), b.getPatternCount());
		for (int i = 0; i < a.getPatternCount(); i++) {
			int [] pattern1 = a.getPattern(i);
			int [] pattern2 = b.getPattern(i);
			for (int j = 0; j < pattern1.length; j++) {
				if (pattern1[j] != pattern2[j]) {
					int h = 3;
					h++;
				}
				assertEquals(pattern1[j], pattern2[j]);
			}
			assertEquals(a.getPatternWeight(i), b.getPatternWeight(i));
		}
	}
}
