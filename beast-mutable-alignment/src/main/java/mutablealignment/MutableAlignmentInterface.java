package mutablealignment;

public interface MutableAlignmentInterface {

	
	public void setSiteValue(int taxonNr, int siteNr, int newValue);
	public int getSiteValue(int taxonNr, int siteNr);


	/**
	 * Set characters in the alignment at a specific taxon for all sites
	 */
	public void setSiteValuesByTaxon(int taxonNr, int [] newValues);
	public int [] getSiteValuesByTaxon(int taxonNr);

	
	/**
	 * Set characters in the alignment at a specific site for all taxa
	 */
	public void setSiteValuesBySite(int siteNr, int [] newValues);
	public int [] getSiteValuesBySite(int siteNr);


	/**
	 * Set all characters in the alignment 
	 */
	public void setSiteValues(int [][] newValues);
	
	/** reset site patterns to previous values **/
	public void resetSitePatterns(int siteNr, int taxonNr, int oldValue);
	public void resetSitePatterns(int siteNr, int[] oldValues);
	public void resetSitePatterns(int[][] oldValues);

}
