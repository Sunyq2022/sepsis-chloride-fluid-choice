# N09 revision code and aggregate results

This package supports the first Critical Care revision dated 24 October 2026. It contains no patient-level files, credentials or source clinical data. Public repository: https://github.com/Sunyq2022/sepsis-chloride-fluid-choice. Initial snapshot tag: revision-2026-10-25.

## Implementations and immutable originals

The archive preserves the original project-relative layout. `02代码/01分析脚本/` contains read-only snapshots of the retained implementation and dependencies. New revision code is under `03结果/21_FULL_ICU_CAUSAL_PROGRAM_20260907/N09/critical_care_revision/code/`. Extract the archive into an empty review directory, never over an active research checkout.

The new revision scripts do not call the original `extract()`, `analyse()` or original manuscript builders. They import pure estimator helpers, read retained inputs, and write only to `critical_care_revision`. The original SQL, protocol, manuscript and results were not overwritten. Source checksums are in the integrity JSON files.

## Reproduction levels

1. **Aggregate review without restricted data.** Inspect the supplied numeric CSVs, protocol, source SQL, diagnostic figure source and code. This does not reproduce individual-level estimation.
2. **Re-estimation for an authorized researcher with the frozen private cohort.** Place the authorized `cohort_private.parquet` under `rebuilt_analysis` in the original layout. Use the retained Python environment and run `revision_analysis.py --mode core`. This independently checks all eight main weighted risks, fits the diagnostic model and runs the six-hour landmark with 500 patient bootstraps, seed 20260909. Private cohort data are not distributed.
3. **Targeted raw-source checks.** `revision_analysis.py --mode source` reads MIMIC-IV v3.1 `hosp/labevents.csv.gz` and `icu/inputevents.csv.gz`, using the retained cohort and prior authorized source extracts. Configure the raw-data location in the new script for another host. It checks the 552 trigger records and unspecified prior-fluid items, and writes restricted row caches only under `work/private` (not included here).
4. **Full extraction from raw MIMIC-IV.** The retained extraction SQL/code additionally require `mimiciv_hosp`, `mimiciv_icu`, and the specific MIMIC-derived sepsis3/SOFA/KDIGO/ventilation/vasopressor tables, plus the adult first-ICU source relation. Historical DuckDB cache paths in the original helper no longer exist on this host. The entire original extraction was not rerun in this revision. The archive provides the actual retained extraction implementation, not a claim of a fresh, portable one-command raw-data rebuild. Credentialed readers must configure their own equivalent relations and audit upstream derived-table versions before rebuilding. No alternative definitions should be substituted silently.

## Commands in the existing project environment

```text
python critical_care_revision/code/revision_analysis.py --mode core
python critical_care_revision/code/revision_analysis.py --mode source
python critical_care_revision/code/revision_figures.py
```

Run from the N09 context using the original relative directory layout. DOCX authoring uses the separately bundled document runtime. `render_revision.ps1` exports only the new revision documents to internal QA PDFs/PNGs and closes only those read-only documents, not other user documents. Original rendering through the packaged LibreOffice helper was unavailable on this Windows host; WPS export was used and is documented.

## Interpretation and audit boundary

The original main contrast is LR minus saline. Revision tables reorient it to saline minus LR and saline divided by LR; confidence limits are reversed when reciprocated. The no-identified-prior-LR analysis means the 24 hours **before the chloride trigger**, not lifelong LR-naive status. One patient received identified LR between the trigger and time zero; other fluids had unspecified composition.

The landmark selects patients alive, in the index ICU and composite-event-free at six hours, then follows to the original day 7. It changes the target population and conditions on post-treatment status. It is not a correction proving causality. The composite interval includes no difference.

The script's runnable checks validate original-risk reproduction, fixed cohort size, temporal source consistency, 500 successful landmark bootstraps, and original-file integrity. They do not validate clinical exchangeability or independent KRT adjudication.

## Versioning and corrections

This repository is a versioned research record, not a claim that the analyses are error-free. The initial public tag preserves the results available on 25 October 2026. Corrections will be committed on main with a CHANGELOG entry describing the reason, validation, and affected estimates or interpretation. Published tags will not be moved. Cite a commit or tag for each manuscript version. If an error changes a submitted or published conclusion, the manuscript and author response must also be corrected. Previous positive results are not guaranteed to remain unchanged.

## Licensing

No new blanket reuse license is assigned to project-authored code in this deposit. Public availability alone does not grant unrestricted reuse. MIMIC-derived SQL retains its upstream licensing; see THIRD_PARTY_NOTICES.md. Data licensing and credentialed-access conditions are separate; no patient-level or source clinical data are distributed.
