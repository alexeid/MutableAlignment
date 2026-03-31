# MutableAlignment

Alignments that can change during MCMC for [BEAST 3](https://github.com/CompEvol/beast3).

## Maven coordinates

```xml
<dependency>
    <groupId>io.github.compevol</groupId>
    <artifactId>mutable-alignment</artifactId>
    <version>0.1.0-SNAPSHOT</version>
</dependency>
```

JPMS module: `mutable.alignment`

## API

Replace the standard `Alignment` with `mutablealignment.MutableAlignment` and `TreeLikelihood` with `mutablealignment.MATreeLikelihood`.

The `MutableAlignment` can be changed by operators using:

* `setSiteValue()` — set a single character in the alignment
* `setSiteValuesByTaxon()` — set a sequence for a taxon
* `setSiteValuesBySite()` — set site values for all taxa (single site)
* `setSiteValues()` — set the whole alignment

Note that `setSiteValue()` and `setSiteValuesByTaxon()` are inefficiently implemented: partials will be recalculated for all sites. A future BEAGLE-aware implementation will probably suffer from the same limitation due to the BEAGLE API not allowing efficient updates of single sites.
