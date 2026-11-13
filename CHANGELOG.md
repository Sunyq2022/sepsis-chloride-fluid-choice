# Second-round revision (2026-11-13)

- Excluded two pre-bolus creatinine stage 3 cases and corrected version-specific dialysis codes.
- Refit eight comparisons with 500 patient-level resamples each.
- Added source relation SQL and 31 derived definitions with content hashes.
- Reproduced the corrected comparison directly from raw creatinine records with the retained 552-person input, without the historical case-review caches.
- The original public tag remains fixed.

# Change log

## 2026-10-25 — Participant flow documentation

Added a GitHub-rendered participant flowchart, sequential exclusion table, source links and six-hour landmark selection accounting. Existing code, estimates and the revision-2026-10-25 tag are unchanged.


## revision-2026-10-25

Initial public deposition of retained extraction/analysis code and aggregate results. No analysis rerun for deposition. Overall cohort: 552 (371 NS, 181 LR). Saline minus LR kidney-composite risk difference: 6.35 percentage points (95% CI 1.10 to 12.08). No-prior-LR and six-hour-landmark intervals include zero. These associations do not establish causal protection.

Removed two upstream SQL example comments containing record identifiers from public copies only; executable SQL unchanged. Historical Code_README.md is retained for provenance; root README supplies current deposition status.

Future corrections must document reasons, checks and affected findings; existing tagged snapshots remain available.
