"""Build all adult first-ICU source concepts; no causal effect is estimated here.

Reuses official local MIMIC concepts, with dependency closure and resumable tables.
Raw ICU/lab records retain charttime/storetime; missing phenotype is not no AKI.
"""
from pathlib import Path
import argparse
import csv
import gzip
import hashlib
import json
import re
import time
import duckdb

ROOT = Path(__file__).resolve().parents[2]
RAW = Path('E:/ICUdatabase/MIMIC/mimic-iv')
CONCEPT = Path('E:/ICUdatabase/MIMIC/mimic-code-main/mimic-iv/concepts_duckdb')
OUT = ROOT / '03结果/14_SA_AKI_HYDROCORTISONE_20260906'
TARGETS = ['sepsis3', 'kdigo_stages', 'rrt', 'charlson', 'vitalsign',
           'chemistry', 'complete_blood_count', 'bg', 'weight_durations',
           'norepinephrine', 'vasopressin', 'ventilation', 'icustay_detail']

def sqlstr(p):
    return str(p).replace('\\', '/').replace("'", "''")

def dependencies():
    files = {}
    for p in sorted(CONCEPT.rglob('*.sql')):
        s = p.read_text(encoding='utf-8')
        m = re.search(r'CREATE TABLE mimiciv_derived\.(\w+)', s, re.I)
        if m:
            files[m[1]] = (p, s)
    ordered, seen = [], set()
    def visit(name):
        if name in seen:
            return
        seen.add(name)
        p, s = files[name]
        for dep in sorted(set(re.findall(r'mimiciv_derived\.(\w+)', s)) - {name}):
            visit(dep)
        ordered.append((name, p, s))
    for name in TARGETS:
        visit(name)
    return ordered

def concepts_to_rebuild(ordered, previous):
    """Invalidate downstream materialized concepts when an upstream SQL changes."""
    changed = set()
    for name, _, sql in ordered:
        parents = set(re.findall(r'mimiciv_derived\.(\w+)', sql, re.I)) - {name}
        digest = hashlib.sha256(sql.encode()).hexdigest()
        if (previous.get(name, {}).get('sha256') != digest
                or previous.get(name, {}).get('invalidated', False)
                or parents & changed):
            changed.add(name)
    return changed

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--check', action='store_true')
    a = ap.parse_args()
    ordered = dependencies()
    raw_tables = sorted(set(re.findall(r'mimiciv_(hosp|icu)\.(\w+)', '\n'.join(s for _, _, s in ordered))))
    if a.check:
        assert len({n for n, _, _ in ordered}) == len(ordered)
        assert all((RAW / schema / (table + '.csv.gz')).is_file() for schema, table in raw_tables)
        print(json.dumps({'concepts': len(ordered), 'source_tables': raw_tables, 'check': 'PASS'}))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    ledger_path = OUT / 'source_build.json'
    ledger = json.loads(ledger_path.read_text()) if ledger_path.exists() else {'status': 'RUNNING', 'tables': {}, 'concepts': {}, 'effect_estimated': False}
    con = duckdb.connect(str(OUT / 'sa_aki.duckdb'))
    con.execute("SET memory_limit='4GB'")
    con.execute('SET threads=2')
    con.execute('SET preserve_insertion_order=false')
    con.execute(f"SET temp_directory='{sqlstr(OUT / 'duckdb_temp')}'")
    def save():
        ledger_path.write_text(json.dumps(ledger, indent=2), encoding='utf-8')
    def run_table(schema, table, where=''):
        key = f'{schema}.{table}'
        if key in ledger['tables']:
            return
        start = time.monotonic()
        print('READ', key, flush=True)
        p = RAW / schema / (table + '.csv.gz')
        con.execute(f"CREATE OR REPLACE TABLE mimiciv_{schema}.{table} AS SELECT * FROM read_csv_auto('{sqlstr(p)}'){where}")
        n = con.execute(f'SELECT count(*) FROM mimiciv_{schema}.{table}').fetchone()[0]
        ledger['tables'][key] = {'rows': n, 'source_bytes': p.stat().st_size, 'source_mtime_ns': p.stat().st_mtime_ns, 'seconds': round(time.monotonic()-start, 2), 'filter': where}
        con.execute('CHECKPOINT')
        save()
        print('DONE', key, n, flush=True)
    try:
        for schema in ['hosp', 'icu', 'derived']:
            con.execute(f'CREATE SCHEMA IF NOT EXISTS mimiciv_{schema}')
        con.execute('CREATE SCHEMA IF NOT EXISTS study')
        for schema, table in [('hosp','patients'),('hosp','admissions'),('icu','icustays')]:
            run_table(schema, table)
        con.execute('''CREATE OR REPLACE TABLE study.first_icu AS
          WITH ranked AS (SELECT i.*,p.gender,p.anchor_age + year(i.intime)-p.anchor_year AS age,
          row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS rn
          FROM mimiciv_icu.icustays i JOIN mimiciv_hosp.patients p USING(subject_id))
          SELECT * FROM ranked WHERE rn=1 AND age>=18''')
        ledger['adult_first_icu_n'] = con.execute('SELECT count(*) FROM study.first_icu').fetchone()[0]
        save()
        # Filter only item-bearing large tables. All referenced item constants are retained.
        itemids = sorted(set(re.findall(r'\b[25]\d{4,5}\b', '\n'.join(s for _,_,s in ordered))))
        for schema, table in raw_tables:
            if f'{schema}.{table}' in ledger['tables']:
                continue
            with gzip.open(RAW/schema/(table+'.csv.gz'), 'rt', encoding='utf-8') as f:
                header = next(csv.reader(f))
            key = next((k for k in ['stay_id','hadm_id','subject_id'] if k in header), None)
            where = f' WHERE {key} IN (SELECT {key} FROM study.first_icu)' if key else ''
            if table=='labevents':
                # Include preadmission results with null hadm_id; do not lose renal history.
                where = ' WHERE subject_id IN (SELECT subject_id FROM study.first_icu)'
            if table in ['chartevents','labevents']:
                where += ' AND itemid IN ('+','.join(itemids)+')'
            run_table(schema, table, where)
        rebuild = concepts_to_rebuild(ordered, ledger['concepts'])
        # Persist descendants before executing SQL so a failed run cannot resume
        # with rebuilt parents but stale, apparently unchanged children.
        for name in rebuild:
            ledger['concepts'].setdefault(name, {})['invalidated'] = True
        save()
        for name, p, s in ordered:
            digest = hashlib.sha256(s.encode()).hexdigest()
            if name not in rebuild:
                continue
            print('CONCEPT', name, flush=True)
            start = time.monotonic()
            con.execute(s)
            ledger['concepts'][name] = {'sha256':digest,'file':str(p),'rows':con.execute(f'SELECT count(*) FROM mimiciv_derived.{name}').fetchone()[0], 'seconds':round(time.monotonic()-start,2)}
            con.execute('CHECKPOINT')
            save()
        ledger['status']='SOURCE_CONCEPTS_BUILT_NOT_ANALYSIS_COHORT'
        save()
        print(ledger['status'], flush=True)
    except Exception as e:
        ledger['status']='FAILED';ledger['error']=repr(e);save();raise
    finally:
        con.close()

if __name__=='__main__':
    main()
