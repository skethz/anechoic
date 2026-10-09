# Addendum 4.1: the SB methods on TSP with their papers' treatment of local fields (8 October 2026, 14:42 CEST)

## What changes

- **Before.** PROTOCOL.md gave every method the field b_i the same way. For the SB methods this was b added to every J·x product, i.e. a frozen ancillary spin x_0 = 1 (our choice).
- **What the paper prescribes.** Goto et al. 2021 (main text, and its section S1) handle local fields differently: "introducing an ancillary spin reduces the Ising problem to the one without local fields". The ancillary spin is a dynamical spin of the field-free (N+1)-spin problem, and the solution maps back by s_i·s_N.
- **Change.** Under the rule "strictly follow what they describe", aSB, bSB and dSB on TSP (the only problem with fields) are re-run with the ancillary-spin reduction (`ancilla.py`).
  - The normalisations then follow from the (N+1)-spin matrix J': c0 = 0.5/(⟨J'⟩√(N+1)) and ξ0 = 0.7/(SD(J')√(N+1)).
  - bSB and dSB Δt is chosen by the paper's procedure.
  - aSB is run at Δt = 0.9 (the paper's K2000 value) and Δt = 0.5 (our choice, as in Addendum 4).
- **Check.** The reduction is checked exactly on all six TSP instances: E'(s, +1) = E'(−s, −1) = E(s), and decoding returns s (`run_a4b.py check`).
- **Unaffected:**
  - **No fields:** Max-Cut and GPP.
  - **Native fields:** Neal, STATICA, APC-SCA and ReAIM all include h in their local fields or Hamiltonians.
  - **TEC:** its paper does not treat fields; on TSP it is under "settings not specified".

## Seeds and files (fresh)

- **Finals:** SeedSequence([20261008, 54, p, i, m, s]) for bSB/dSB; 56 for aSB at Δt 0.9; 57 for aSB at Δt 0.5.
- **bSB/dSB Δt pilots:** [.., 55, .., k].
- **Output:** `results_a4/final/sb_anc/`.
- **Earlier results.** The frozen-x_0 SB results on TSP (A1, A3) are kept, labelled "frozen ancilla (our earlier choice)".

## Disclosure

- **Timing.** Written at 14:42 CEST, while the Addendum 4 validation runs were in progress. The Addendum 4 Table 8b re-runs had finished: their TEC/aSB summaries, including aSB diverging on TSP, had been seen.
- **What prompted it.** Re-reading the 2021 paper's treatment of fields while diagnosing aSB on TSP.
