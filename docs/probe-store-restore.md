# Probe APIs and the store/restore lifecycle

In this document, **probe** means a speculative call to
`MATreeLikelihood.getLogProbsForStateSequence(nodeNr, sites)` (or its
partials variant `getLogProbsForPartialsSequence`): the caller is asking
*what would the per-site log likelihoods be if leaf `nodeNr` had this
sequence?* without committing to it as the new state. Gibbs-style
operators issue many such probes per proposal to build up their proposal
distribution.

The probe walks from `nodeNr` up to the root, recomputing partials at
every ancestor — which writes into BEAST2's `LikelihoodCore` buffers.
That side effect is the source of all the subtleties below: it has to be
managed against BEAST2's store/restore lifecycle so that the cached
partials representing the current accepted state survive untouched.

## The buffer model

`BeerLikelihoodCore` double-buffers internal-node partials with two
per-node arrays and a `currentPartialsIndex[a]` (call it `c`) that selects
the live slot, plus a `storedPartialsIndex[a]` (call it `s`) used by
`restore()`. After `store()` (which copies `c → s`), both `c` and `s` point
to the same data slot — that's the "stored" state to roll back to.

`setNodePartialsForUpdate(a)` flips `c[a]`. The contract is: callers flip
*before* writing new partials, so writes land in the scratch slot and the
stored slot stays intact. `restore()` swaps the index arrays, which has
the effect of pointing `c` back at the stored slot.

Tip states (when `useAmbiguities=false`) are *not* double-buffered —
there's a single `int[] states[nodeIndex]`. `MATreeLikelihood` patches over
this with `tempTipNodes` (track which tips were probed) and resync from
`MutableAlignment` in `restore()`. This piece is unavoidable.

## The MCMC iteration timing

In one MCMC iteration, the framework calls (in order):

1. `state.store(sample)` — bookkeeping only, doesn't touch CalculationNodes
2. `operator.proposal()` — *probes happen here*
3. `state.storeCalculationNodes()` → `MATreeLikelihood.store()` — **after** the proposal, **before** evaluation
4. `posterior.calculateLogP()` → `MATreeLikelihood.calculateLogP()`
5. accept → `acceptCalculationNodes()`, or reject → `state.restore()` + `restoreCalculationNodes()`

Step 3 is the trap of an earlier design (see history below). Anything
`store()` clears is wiped *between probes and calculateLogP* — so any
probe-tracking state set up during step 2 must either (a) survive `store()`
or (b) be eliminated entirely before `store()` is called. The current code
takes path (b) for partials.

## Current design (commit `66e8b99` and later): hermetic probes

`calcPatternLogLikelihoods` is now hermetic: it flips each touched
ancestor's partials index to the scratch slot before writing, computes
pattern log likelihoods at the root, then flips each ancestor back. After
the probe call returns, the partials indices are exactly what they were
on entry — the probe leaves *no trace* in the `LikelihoodCore`'s buffer
indices for the framework's store/restore to mishandle.

```mermaid
sequenceDiagram
    participant MCMC
    participant Op as Operator
    participant MA as MATreeLikelihood
    participant LC as likelihoodCore

    Note over LC: partials[0]=A (current), partials[1]=junk<br/>c=0, s=0

    MCMC->>Op: proposal()
    Op->>MA: getLogProbsForStateSequence(...)
    MA->>LC: setNodePartialsForUpdate(a) [c: 0→1]
    MA->>LC: calculatePartials → writes partials[1][a]
    Note over LC: partials[0]=A, partials[1]=probe<br/>c=1, s=0
    MA->>LC: setNodePartialsForUpdate(a) [c: 1→0]
    Note over LC: partials[0]=A, partials[1]=probe orphaned<br/>c=0, s=0 (back to entry)

    MCMC->>MA: store ok, no partials state to preserve

    MCMC->>MA: calculateLogP()
    MA->>LC: super.calculateLogP → traverse flips c again [c: 0→1]<br/>writes new state B into partials[1]
    Note over LC: partials[0]=A (intact!), partials[1]=B<br/>c=1, s=0

    alt accept
        MCMC->>MA: accept → LC stores: s := c
        Note over LC: partials[1]=B, c=s=1
    else reject
        MCMC->>LC: restore swaps c and s
        Note over LC: c=0, s=1 -- reads partials[0]=A
    end
```

For tip states, `tempTipNodes` is still needed and `store()` must not
clear it — that is the one cross-method invariant left.

## Earlier design (commits `0a33bb7`/`c99c8f8`/`af09f09`)

Earlier versions tracked flipped ancestors in a class-level
`flippedNodes` Set and undid the flips in an explicit loop at the start
of `calculateLogP()`. That worked in principle, but commits `0a33bb7` and
`c99c8f8` cleared `flippedNodes` inside `store()` — wiping the tracking
between probes and the undo loop. The framework's traverse then flipped
`c` a *second* time, rotating it back onto the stored slot and
overwriting the genuine pre-store partials. The bug surfaced as a
likelihood mismatch in `phylonco`'s `ExchangeGibbsOperator` (cached
−541.026 vs fresh −570.892, Δ 29.87).

```mermaid
sequenceDiagram
    participant MCMC
    participant Op as Operator
    participant MA as MATreeLikelihood
    participant LC as likelihoodCore

    Note over LC: partials[0]=A (current), partials[1]=junk<br/>c=0, s=0

    MCMC->>Op: proposal()
    Op->>MA: getLogProbsForStateSequence(...)
    MA->>LC: setNodePartialsForUpdate(a) [c: 0→1]
    Note right of MA: flippedNodes += a
    MA->>LC: calculatePartials → writes partials[1][a]
    Note over LC: partials[0]=A, partials[1]=probe<br/>c=1, s=0

    rect rgb(255, 200, 200)
    MCMC->>MA: store clears flippedNodes (BUG)
    end

    MCMC->>MA: calculateLogP()
    Note right of MA: undo loop is empty -- c stays at 1
    MA->>LC: super.calculateLogP → traverse flips c again [c: 1→0]<br/>writes new state B into partials[0]
    Note over LC: partials[0]=B (overwrote A!), partials[1]=probe<br/>c=0, s=0

    MCMC->>LC: restore swaps c and s
    Note over LC: partials[0] holds the corrupt value B instead of A
```

`af09f09` patched the immediate symptom by removing the two `clear()`
calls in `store()`. `66e8b99` then refactored to the hermetic design
above, eliminating the cross-method partials state altogether.

## Why this is the cleanest available design

- **Tip states must use `tempTipNodes`.** `setNodeStatesForUpdate` is a
  no-op in `BeerLikelihoodCore`, so probes that mutate leaf states have
  to resync from `MutableAlignment` in `restore()` regardless of what
  partials handling we choose. BEAGLE's tip states aren't double-buffered
  either, so the `BeagleMATreeLikelihood` path needs the same workaround.
- **Partials don't need cross-method state.** The "flip up, compute,
  flip back" pattern is local to `calcPatternLogLikelihoods` and leaves
  no cross-method invariants (other than the one tip-state set we already
  need). No `store()`/`accept()`/`restore()` machinery to keep in sync;
  no don't-clear-in-store trap to remember.
- **BEAGLE's rescale recursion is handled** by threading the per-call
  `flipped` set through the recursive call, so each ancestor is still
  flipped at most once across rescale attempts and flipped back exactly
  once on the final success path.

## Edge case left open

If an operator uses probes and then returns `Double.NEGATIVE_INFINITY`,
MCMC's failure path calls `state.restore()` but skips
`restoreCalculationNodes()` (when the operator's
`requiresStateInitialisation()` is `true`, the default). `MA.restore()`
never runs, so `tempTipNodes` is not cleared and the tip states are not
resynced.

For partials this is now harmless — the hermetic flip-up/flip-back
already happened inside the probe call, so the buffer indices are clean.
For tip states, the next probe in the next iteration would write fresh
states (overwriting the stale ones) before computing, so single-probe
operators still get the right answer. Cross-iteration state leakage of
`tempTipNodes` is a small leak, not a correctness bug.

This isn't triggered by `phylonco`'s `ExchangeGibbsOperator` (all its
`NEGATIVE_INFINITY` early-returns happen *before* the first probe) but
worth knowing if any future probe-using operator returns
`NEGATIVE_INFINITY` post-probe.
