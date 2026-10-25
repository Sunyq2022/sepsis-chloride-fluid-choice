SELECT l.labevent_id,l.itemid,l.charttime,l.storetime,l.valuenum FROM read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/hosp/labevents.csv.gz') l SEMI JOIN cohort c USING(labevent_id);

SELECT c.stay_id,c.exposure,c.lr_ml,c.trigger_time,c.t0,i.itemid,i.starttime,i.endtime,i.amount,i.amountuom,i.ordercategoryname
    FROM read_csv_auto('E:/ICUdatabase/MIMIC/mimic-iv/icu/inputevents.csv.gz') i JOIN cohort c USING(stay_id)
    WHERE i.itemid IN (226361,226364,226375) AND i.starttime<c.t0
    AND i.statusdescription IS DISTINCT FROM 'Rewritten';
