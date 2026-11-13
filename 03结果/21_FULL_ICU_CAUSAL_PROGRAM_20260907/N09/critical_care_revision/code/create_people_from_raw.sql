-- MIMIC-IV 3.1; execute in the source database attached as s by the analysis.
CREATE OR REPLACE TABLE people AS
WITH ranked AS (
 SELECT i.*,p.gender,p.anchor_age+year(i.intime)-p.anchor_year AS age,
 row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS rn
 FROM read_csv_auto('E:\ICUdatabase\MIMIC\mimic-iv/icu/icustays.csv.gz') i
 JOIN read_csv_auto('E:\ICUdatabase\MIMIC\mimic-iv/hosp/patients.csv.gz') p USING(subject_id)
), first_icu AS (SELECT * FROM ranked WHERE rn=1 AND age>=18)
SELECT f.*,a.admittime,a.dischtime,a.deathtime,a.admission_type,a.hospital_expire_flag,p.anchor_year_group
FROM first_icu f JOIN read_csv_auto('E:\ICUdatabase\MIMIC\mimic-iv/hosp/admissions.csv.gz') a ON f.hadm_id=a.hadm_id
JOIN read_csv_auto('E:\ICUdatabase\MIMIC\mimic-iv/hosp/patients.csv.gz') p ON f.subject_id=p.subject_id;
