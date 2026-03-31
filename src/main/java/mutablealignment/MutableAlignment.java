package mutablealignment;



import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.w3c.dom.Node;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.datatype.DataType;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;

@Description("Alignment that can be sampled by MCMC")
public class MutableAlignment extends Alignment implements MutableAlignmentInterface {
	
	public MutableAlignment() {
	}
	
	public MutableAlignment(Alignment other) {
		for (Input<?> input : other.listInputs()) {
			setInputValue(input.getName(), input.get());
		}
		initAndValidate();
	}
	
	
	@Override
	public void initAndValidate() {
		// sort sequences alphabetically
		// TODO why are we doing this? What if the names are integers?
		//Collections.sort(sequenceInput.get(), (o1,o2) -> {
		//	return o1.taxonInput.get().compareTo(o2.taxonInput.get());
		//});
		
		super.initAndValidate();
	}
	
	/**
	 * Set 1 character in the alignment at a specific taxon and site
	 */
	public void setSiteValue(int taxonNr, int siteNr, int newValue) {
		startEditing(null);
		editList.add(new Edit(siteNr, taxonNr, sitePatterns[siteNr][taxonNr], newValue));
		sitePatterns[siteNr][taxonNr] = newValue;
	}

	public int getSiteValue(int taxonNr, int siteNr) {
		return sitePatterns[siteNr][taxonNr];
	}


	/**
	 * Set characters in the alignment at a specific taxon for all sites
	 */
	public void setSiteValuesByTaxon(int taxonNr, int [] newValues) {
		startEditing(null);
		int [] oldValues = new int[newValues.length];
		for (int i = 0; i < newValues.length; i++) {
			oldValues[i] = sitePatterns[i][taxonNr]; 
		}
		editList.add(new Edit(oldValues, taxonNr, newValues));
		for (int i = 0; i < newValues.length; i++) {
			sitePatterns[i][taxonNr] = newValues[i];
		}
	}

	public int [] getSiteValuesByTaxon(int taxonNr) {
		int [] seq = new int[sitePatterns.length];
		for (int i = 0; i < seq.length; i++) {
			seq[i] = sitePatterns[i][taxonNr];
		}
		return seq;
	}

	@Override
	public void startEditing(Operator operator) {
		super.startEditing(operator);
		setSomethingIsDirty(true);
	}
	
	/**
	 * Set characters in the alignment at a specific site for all taxa
	 */
	public void setSiteValuesBySite(int siteNr, int [] newValues) {
		startEditing(null);
		editList.add(new Edit(siteNr, sitePatterns[siteNr].clone(), newValues));
		System.arraycopy(newValues, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
	}

	public int [] getSiteValuesBySite(int siteNr) {
		return getPattern(siteNr);
	}


	/**
	 * Set all characters in the alignment 
	 */
	public void setSiteValues(int [][] newValues) {
		startEditing(null);
		int n = newValues[0].length;
		int [][] oldValues = new int[newValues.length][n];
		for (int i = 0; i < oldValues.length; i++) {
			System.arraycopy(sitePatterns[i], 0, oldValues[i], 0, n);
		}
		editList.add(new Edit(oldValues, newValues));
		for (int i = 0; i < newValues.length; i++) {
			System.arraycopy(newValues[i], 0, sitePatterns[i], 0, sitePatterns[i].length);
		}
	}
	
	
	
	/** reset site patterns to previous values **/
	public void resetSitePatterns(int siteNr, int taxonNr, int oldValue) {
		sitePatterns[siteNr][taxonNr] = oldValue;
	}

	public void resetSitePatterns(int siteNr, int[] oldValues) {
		System.arraycopy(oldValues, 0, sitePatterns[siteNr], 0, sitePatterns[siteNr].length);
	}

	public void resetSitePatterns(int[][] oldValues) {
		for (int i = 0; i < oldValues.length; i++) {
			System.arraycopy(oldValues[i], 0, sitePatterns[i], 0, sitePatterns[i].length);
		}
	}
	
	
	
	/** Loggable Interface **/
	@Override
	public void init(PrintStream out) {
		for (int i = 0; i < taxaNames.size(); i++) {
			out.print(taxaNames.get(i) + "\t");
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		for (int i = 0; i < taxaNames.size(); i++) {
			for (int j = 0; j < sitePatterns.length; j++) {
				out.print(getDataType().getCharacter(sitePatterns[j][i]));
			}
			out.print("\t");
		}
	}
	
	@Override
	public void close(PrintStream out) {
	}
	

	/** StateNode stuff **/
	
	enum EditType {singleSite, allTaxa, allSites, all}
	
	protected List<Edit> editList = new ArrayList<>();
	
	@Override
	protected void store() {
		editList.clear();
		hasStartedEditing = false;
		super.store();
	}

	
	@Override
	protected void accept() {
		hasStartedEditing = false;
		editList.clear();
		super.accept();
	}
	
	@Override
	public void restore() {
		
		for (int i = editList.size()-1; i>=0; i--) {
			editList.get(i).undo(this);
		}
		editList.clear();
		hasStartedEditing = false;
		super.restore();
	}

	/**
	 * calculate patterns from sequence data *
	 */
	protected void calcPatterns(boolean log) {
		int taxonCount = counts.size();
		int siteCount = counts.get(0).size();

		// convert data to transposed int array
		sitePatterns = new int[siteCount][taxonCount];
		for (int i = 0; i < taxonCount; i++) {
			List<Integer> sites = counts.get(i);
			for (int j = 0; j < siteCount; j++) {
				sitePatterns[j][i] = sites.get(j);
			}
		}
		int patterns = siteCount;
		
		// reserve memory for patterns
		patternWeight = new int[patterns];
		Arrays.fill(patternWeight, 1);

		// find patterns for the sites
		patternIndex = new int[siteCount];
		for (int i = 0; i < siteCount; i++) {
			patternIndex[i] = i;
		}

		if (siteWeights != null) {
			Arrays.fill(patternWeight, 0);
			for (int i = 0; i < siteCount; i++) {
				patternWeight[patternIndex[i]] += siteWeights[i];
			}
		}

		// determine maximum state count
		// Usually, the state count is equal for all sites,
		// though for SnAP analysis, this is typically not the case.
		maxStateCount = 0;
		for (int m_nStateCount1 : stateCounts) {
			maxStateCount = Math.max(maxStateCount, m_nStateCount1);
		}
		// report some statistics
		if (log && taxaNames.size() < 30) {
			for (int i = 0; i < taxaNames.size(); i++) {
				Log.info.println(taxaNames.get(i) + ": " + counts.get(i).size() + " " + stateCounts.get(i));
			}
		}

		if (stripInvariantSitesInput.get()) {
			// don't add patterns that are invariant, e.g. all gaps
			if (log)
				Log.info.println("Stripping invariant sites");

			int removedSites = 0;
			for (int i = 0; i < patterns; i++) {
				int[] pattern = sitePatterns[i];
				int value = pattern[0];
				boolean isInvariant = true;
				for (int k = 1; k < pattern.length; k++) {
					if (pattern[k] != value) {
						isInvariant = false;
						break;
					}
				}
				if (isInvariant) {
					removedSites += patternWeight[i];
					patternWeight[i] = 0;

					if (log)
						Log.info.print(" <" + value + "> ");
				}
			}
			if (log)
				Log.info.println(" removed " + removedSites + " sites ");
		}
	} // calcPatterns

	public Integer [] getDirtySequenceIndices() {
		Set<Integer> dirtySequences = new HashSet<>();
		for (Edit edit : editList) {
			switch (edit.type) {
			case all:
			case allTaxa:
				for (int i = 0; i < getTaxonCount(); i++) {
					dirtySequences.add(i);
				}
				break;
			case allSites:
			case singleSite:
				dirtySequences.add(edit.taxonNr);
				break;
			}
		}
		return dirtySequences.toArray(new Integer[] {});
	}

	@Override
	/* 
	 * be aware the toString() is used to restore from the state,
	 * so if you change this, then fromXML() needs to be updated as well
	 */
	public String toString() {
		int siteCount = getSiteCount();
		int taxonCount = getTaxonCount();
		
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < taxonCount; i++) {
			b.append(getTaxaNames().get(i) + ": ");
			for (int j = 0; j < siteCount; j++) {
				b.append(getDataType().getCharacter(getPattern(i, j)));
			}
			b.append("\n");
		}
		return b.toString();
	}

	@Override
	public String toXML() {
        return "<statenode id='" + normalise(getID()) + "'>" +
                normalise(toString()) +
                "</statenode>\n";
    }

    /** ensure XML identifiers get proper escape sequences **/
    private String normalise(String str) {
    	if (str == null) {
    		return null;
    	}
    	str = str.replaceAll("&", "&amp;");    	
    	str = str.replaceAll("'", "&apos;");
    	str = str.replaceAll("\"", "&quot;");
    	str = str.replaceAll("<", "&lt;");
    	str = str.replaceAll(">", "&gt;");
    	return str;
    }

	@Override
	public void fromXML(Node node) {
		String str = node.getTextContent();
		String [] strs = str.split("\n");
		for (String s : strs) {
			String []strs2 = s.split(":");
			String taxon = strs2[0];
			String sequences = strs2[1].strip();
			DataType dataType = getDataType();
			int [] seq = dataType.stringToEncoding(sequences).stream().mapToInt(i->i).toArray();
			int taxonNr = getTaxonIndex(taxon);
			for (int siteNr = 0; siteNr < seq.length; siteNr++) {
				sitePatterns[siteNr][taxonNr] = seq[siteNr];
			}
		}
	}

	@Override
	public StateNode copy() {
		return new MutableAlignment(this);
	}
	
	@Override
	public void assignFromFragile(StateNode other) {
		MutableAlignment src = (MutableAlignment)other;
		if (this.getPatternCount() != src.getPatternCount()) {
			throw new IllegalArgumentException("assignFromFragile() Expected replacement to be of equal number of patterns (this is realy fragile)");
		}
		for (int i = 0; i < src.getPatternCount(); i++) {
			resetSitePatterns(i, src.getPattern(i));			
		}
	}

}
