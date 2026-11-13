"""Reconstruct the missing people relation and check retained source SQL hashes."""
from pathlib import Path
import json,hashlib,shutil
import duckdb,pandas as pd
import revision_analysis as a
R=a.REV;O=R/'results/round2';RAW=Path('E:/ICUdatabase/MIMIC/mimic-iv')
sql=f"""-- MIMIC-IV 3.1; execute in the source database attached as s by the analysis.
CREATE OR REPLACE TABLE people AS
WITH ranked AS (
 SELECT i.*,p.gender,p.anchor_age+year(i.intime)-p.anchor_year AS age,
 row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS rn
 FROM read_csv_auto('{RAW}/icu/icustays.csv.gz') i
 JOIN read_csv_auto('{RAW}/hosp/patients.csv.gz') p USING(subject_id)
), first_icu AS (SELECT * FROM ranked WHERE rn=1 AND age>=18)
SELECT f.*,a.admittime,a.dischtime,a.deathtime,a.admission_type,a.hospital_expire_flag,p.anchor_year_group
FROM first_icu f JOIN read_csv_auto('{RAW}/hosp/admissions.csv.gz') a ON f.hadm_id=a.hadm_id
JOIN read_csv_auto('{RAW}/hosp/patients.csv.gz') p ON f.subject_id=p.subject_id;
"""
(R/'code/create_people_from_raw.sql').write_text(sql,encoding='utf-8')
c=duckdb.connect();c.execute(sql)
d=pd.read_parquet(a.OLD/'cohort_private.parquet');c.register('retained',d)
assert c.execute('select count(*) from people').fetchone()[0]==65366
assert c.execute('select count(*) from retained r join people p using(subject_id,hadm_id,stay_id) where r.age=p.age and r.intime=p.intime').fetchone()[0]==552
ledger=json.loads((a.ROOT/'03结果/14_SA_AKI_HYDROCORTISONE_20260906/source_build.json').read_text('utf-8'))
concept=Path('E:/ICUdatabase/MIMIC/mimic-code-main/mimic-iv/concepts_duckdb')
dest=R/'code/derived_sql';dest.mkdir(exist_ok=True)
rows=[]
for name,entry in ledger['concepts'].items():
    matches=list(concept.rglob(name+'.sql'))
    if len(matches)!=1:raise RuntimeError(name)
    p=matches[0];digest=hashlib.sha256(p.read_text('utf-8').encode()).hexdigest()
    rows.append({'table':name,'path':str(p.relative_to(concept)),'build_sha256':entry.get('sha256'),'current_sha256':digest,'matches_build':entry.get('sha256')==digest})
    if entry.get('sha256')==digest:
        target=dest/p.relative_to(concept);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
pd.DataFrame(rows).to_csv(O/'derived_sql_versions.csv',index=False)
(O/'source_provenance.json').write_text(json.dumps({'mimic_version':'3.1','people_n':65366,'historical_cohort_keys_age_time_matched':552,'derived_tables':len(rows),'all_sql_hashes_match_materialization_ledger':all(x['matches_build'] for x in rows),'upstream_git_commit':'Not recorded; exact SQL SHA-256 identifies retained derived definitions.'},indent=2),encoding='utf-8')
print((O/'source_provenance.json').read_text())
