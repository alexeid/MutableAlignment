package mutablealignment;

import java.io.PrintStream;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.core.Input.Validate;

@Description("Logs whole mutable alignment")
public class MutableAlignmentLogger extends BEASTObject implements Loggable {
	final public Input<MutableAlignment> alignmentInput = new Input<>("alignment", "mutable alignment to be logged", Validate.REQUIRED);

	private MutableAlignment alignment;
	
	@Override
	public void initAndValidate() {
		alignment = alignmentInput.get();
	}
	
	
	@Override
	public void init(PrintStream out) {
		int siteCount = alignment.getSiteCount();
		int taxonCount = alignment.getTaxonCount();
		List<String> taxaNames = alignment.getTaxaNames();
		
		for (int i = 0; i < taxonCount; i++) {
			for (int j = 0; j < siteCount; j++) {
				out.append(taxaNames.get(i) + j + "\t");
			}
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		int siteCount = alignment.getSiteCount();
		int taxonCount = alignment.getTaxonCount();
		
		for (int i = 0; i < taxonCount; i++) {
			for (int j = 0; j < siteCount; j++) {
				out.append(alignment.getPattern(i, j) + "\t");
			}
		}
	}

	@Override
	public void close(PrintStream out) {
	}

}
