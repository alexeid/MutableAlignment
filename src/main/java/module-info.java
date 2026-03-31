open module mutable.alignment {
    requires beast.base;
    requires beast.pkgmgmt;
    requires beagle;
    requires java.xml;

    exports mutablealignment;

    provides beast.base.core.BEASTInterface with
            mutablealignment.MutableAlignment,
            mutablealignment.MATreeLikelihood,
            mutablealignment.BeagleMATreeLikelihood,
            mutablealignment.MutableAlignmentOperator,
            mutablealignment.MutableAlignmentLogger;
}
