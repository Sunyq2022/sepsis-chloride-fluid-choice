# N09 second-round revision code and aggregate results

Second-round Critical Care revision, 13 November 2026. Versioned deposit: https://github.com/Sunyq2022/sepsis-chloride-fluid-choice/tree/revision-2026-11-13. The first-round results remain at the revision-2026-10-25 tag. No patient-level data are distributed.

## Current versus historical files

Current estimates and checks are under critical_care_revision/results/round2. Earlier rebuilt_analysis and critical_care_revision/results files are historical and must not be presented as corrected 550-patient estimates. New code preserves original implementation and frozen endpoints/seed. The clinical eligibility correction excludes two pre-bolus stage 3 specimens returned later; corrected version-specific dialysis codes add no exclusions. Eight comparisons use 500 patient bootstraps each, seed 20260909, with canonical subject/stay ordering.

## Source relations and versions

code/create_people_from_raw.sql creates the missing people relation directly from MIMIC-IV 3.1 patients, admissions and icustays. Set the raw file paths to the credentialed local installation and execute in the source database attached as s downstream. This was verified against 65,366 adult first ICU rows and every historical cohort key, age and ICU entry time.

code/derived_sql contains 31 exact retained SQL definitions. results/round2/derived_sql_versions.csv compares each SHA-256 with the original materialization ledger; all match. An upstream Git commit was not recorded. The manuscript therefore cites content versions rather than an unverified upstream release. Patient-identifier example comments in historical SQL may be removed in the distribution; executable logic is unchanged and the archive manifest records distributed bytes.

## Execution and requirements

The archive retains original project-relative paths. Extract into a separate review directory, never over an active checkout. Raw data require MIMIC credentialing and a valid data-use agreement. The analysis environment has Python 3.13, pandas, NumPy, scikit-learn, DuckDB, pyarrow and threadpoolctl; retained runtime_versions.json documents the original packages.

For authorized re-estimation, first create the people relation in the database attached as s and configure the MIMIC hosp, ICU and retained derived relations in b. Run the retained n09_corrected_estimand.py extract mode in the documented project layout to recreate the historical 552-person cohort and aggregate outcomes. Then run python critical_care_revision/code/round2_analysis.py from the N09 directory. The estimator reads raw MIMIC-IV 3.1 labevents, admissions and diagnoses to determine the additional exclusions and refits all eight comparisons. It requires neither the historical 706-person cohort nor the previous private case-review caches. Adjust the local raw-data path in the script. It writes the corrected private cohort only under work/private/round2. No fresh end-to-end raw event extraction was performed for this deposit.

code/round2_source_provenance.py rebuilds people in memory and validates derived SQL hashes. code/round2_exhibits.py generates the current aggregate diagnostics and plots. Paths are explicitly local and require adaptation elsewhere. Dependencies for full event extraction remain documented in the original source SQL and 31 derived definitions; this revision did not rerun the entire raw-event extraction or claim a portable one-command rebuild. Importing older scripts is not a license to run their entry points.

## Interpretation

The overlap target is not the whole sepsis population. Restricted comparisons target different overlap populations; trigger strata retain residual imbalance. The landmark conditions on post-treatment status and is not a causal correction. Corrected prior-admission codes cannot exclude maintenance dialysis first documented at the index admission. The current aggregate interval excludes zero, but neither this nor individual components establishes a benefit of switching fluids. Historical reconstruction and multiple comparisons remain unaccounted for by ordinary bootstrap intervals.

## Participant flow

See the [current 550-person selection flow](docs/participant_flow_round2.md) and the [historical 552-person flow](docs/participant_flow.md).
