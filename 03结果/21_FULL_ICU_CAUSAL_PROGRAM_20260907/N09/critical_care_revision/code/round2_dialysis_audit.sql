SELECT c.stay_id, dx.hadm_id, dx.icd_code, dx.icd_version,
        dx.hadm_id=c.hadm_id AS index_admission,
        dx.hadm_id<>c.hadm_id AND ad.dischtime<=c.t0 AS completed_prior_admission
        FROM cohort c JOIN read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/admissions.csv.gz') ad ON ad.subject_id=c.subject_id
        JOIN read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/diagnoses_icd.csv.gz') dx ON dx.hadm_id=ad.hadm_id
        WHERE ((dx.icd_version=9 AND replace(upper(dx.icd_code),'.','') IN ('5856','V451','V4511')) OR (dx.icd_version=10 AND replace(upper(dx.icd_code),'.','') IN ('N186','Z992')));
