package mutablealignment;

import mutablealignment.MutableAlignment.EditType;

/** 
 * class for tracking edits to the alignment 
 * used to restore an alignment if necessary
 **/
class Edit {
	EditType type;
	int siteNr;
	int taxonNr;
	Object oldValue;
	Object newValue;
	
	Edit(int siteNr, int taxonNr, int oldValue, int newValue) {
		this.type = EditType.singleSite;
		this.siteNr = siteNr;
		this.taxonNr = taxonNr;
		this.oldValue = oldValue;
		this.newValue = newValue;
	}

	Edit(int siteNr, int [] oldValue, int [] newValue) {
		this.type = EditType.allTaxa;
		this.siteNr = siteNr;
		this.taxonNr = -1;
		this.oldValue = oldValue;
		this.newValue = newValue;
	}

	Edit(int [] oldValue, int taxonNr, int [] newValue) {
		this.type = EditType.allSites;
		this.siteNr = -1;
		this.taxonNr = taxonNr;
		this.oldValue = oldValue;
		this.newValue = newValue;
	}
	
	Edit(int [][] oldValue, int [][] newValue) {
		this.type = EditType.all;
		this.siteNr = -1;
		this.taxonNr = -1;
		this.oldValue = oldValue;
		this.newValue = newValue;
	}

	void undo(MutableAlignment mutableAlignment) {
		switch(type) {
		case singleSite:
			mutableAlignment.resetSitePatterns(siteNr, taxonNr, (int) oldValue);
			break;
		case allTaxa:
			mutableAlignment.resetSitePatterns(siteNr, (int[])oldValue);
			break;
		case allSites:
			int [] oldValues = (int[]) oldValue;
			for (int i = 0; i < oldValues.length; i++) {
				mutableAlignment.resetSitePatterns(i, taxonNr, oldValues[i]);
			}
			break;
		case all:
			mutableAlignment.resetSitePatterns((int[][])oldValue);
			break;
		}
	}
}