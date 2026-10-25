-- Classify oxygen devices and ventilator modes into six clinical categories.

-- Categories include..
--  Invasive oxygen delivery types: 
--      Tracheostomy (with or without positive pressure ventilation) 
--      InvasiveVent (positive pressure ventilation via endotracheal tube,
--          could be oro/nasotracheal or tracheostomy)
--  Non invasive oxygen delivery types (ref doi:10.1001/jama.2020.9524):
--      NonInvasiveVent (non-invasive positive pressure ventilation)
--      HFNC (high flow nasal oxygen / cannula)
--      SupplementalOxygen (all other non-rebreather,
--          facemask, face tent, nasal prongs...)
--  No oxygen device:
--      None

-- When conflicting settings occur (rare), the priority is:
--  trach > mech vent > NIV > high flow > o2

-- Some useful cases for debugging:

-- first we collect all times which have relevant documentation
WITH tm AS (
    SELECT stay_id, charttime
    FROM b.mimiciv_derived.ventilator_setting
    UNION DISTINCT
    SELECT stay_id, charttime
    FROM b.mimiciv_derived.oxygen_delivery
)

, vs AS (
    SELECT tm.stay_id, tm.charttime
        -- source data columns, here for debug
        , o2_delivery_device_1
        , COALESCE(ventilator_mode, ventilator_mode_hamilton) AS vent_mode
        -- case statement determining the type of intervention
        -- done in order of priority: trach > mech vent > NIV > high flow > o2
        , CASE
            -- tracheostomy
            WHEN o2_delivery_device_1 IN
                (
                    'Tracheostomy tube'
                    -- 1135 observations for T-Piece
                    -- could be either InvasiveVent or Tracheostomy, so omit
                    -- 'T-piece',
                    , 'Trach mask ' -- 16435 observations
                )
                THEN 'Tracheostomy'
            -- mechanical / invasive ventilation
            WHEN o2_delivery_device_1 IN
                (
                    'Endotracheal tube'
                )
                OR ventilator_mode IN
                (
                    '(S) CMV'
                    , 'APRV'
                    , 'APRV/Biphasic+ApnPress'
                    , 'APRV/Biphasic+ApnVol'
                    , 'APV (cmv)'
                    , 'Ambient'
                    , 'Apnea Ventilation'
                    , 'CMV'
                    , 'CMV/ASSIST'
                    , 'CMV/ASSIST/AutoFlow'
                    , 'CMV/AutoFlow'
                    , 'CPAP/PPS'
                    , 'CPAP/PSV'
                    , 'CPAP/PSV+Apn TCPL'
                    , 'CPAP/PSV+ApnPres'
                    , 'CPAP/PSV+ApnVol'
                    , 'MMV'
                    , 'MMV/AutoFlow'
                    , 'MMV/PSV'
                    , 'MMV/PSV/AutoFlow'
                    , 'P-CMV'
                    , 'PCV+'
                    , 'PCV+/PSV'
                    , 'PCV+Assist'
                    , 'PRES/AC'
                    , 'PRVC/AC'
                    , 'PRVC/SIMV'
                    , 'PSV/SBT'
                    , 'SIMV'
                    , 'SIMV/AutoFlow'
                    , 'SIMV/PRES'
                    , 'SIMV/PSV'
                    , 'SIMV/PSV/AutoFlow'
                    , 'SIMV/VOL'
                    , 'SYNCHRON MASTER'
                    , 'SYNCHRON SLAVE'
                    , 'VOL/AC'
                )
                OR ventilator_mode_hamilton IN
                (
                    'APRV'
                    , 'APV (cmv)'
                    , 'Ambient'
                    , '(S) CMV'
                    , 'P-CMV'
                    , 'SIMV'
                    , 'APV (simv)'
                    , 'P-SIMV'
                    , 'VS'
                    , 'ASV'
                )
                THEN 'InvasiveVent'
            -- NIV
            WHEN o2_delivery_device_1 IN
                (
                    'Bipap mask ' -- 8997 observations
                    , 'CPAP mask ' -- 5568 observations
                )
                OR ventilator_mode_hamilton IN
                (
                    'DuoPaP'
                    , 'NIV'
                    , 'NIV-ST'
                )
                THEN 'NonInvasiveVent'
            -- high flow nasal cannula
            WHEN o2_delivery_device_1 IN
                (
                    'High flow nasal cannula' -- 925 observations
                )
                THEN 'HFNC'
            -- non rebreather
            WHEN o2_delivery_device_1 IN
                (
                    'Non-rebreather' -- 5182 observations
                    , 'Face tent' -- 24601 observations
                    , 'Aerosol-cool' -- 24560 observations
                    , 'Venti mask ' -- 1947 observations
                    , 'Medium conc mask ' -- 1888 observations
                    , 'Ultrasonic neb' -- 9 observations
                    , 'Vapomist' -- 3 observations
                    , 'Oxymizer' -- 1301 observations
                    , 'High flow neb' -- 10785 observations
                    , 'Nasal cannula'
                )
                THEN 'SupplementalOxygen'
            WHEN o2_delivery_device_1 IN
                (
                    'None'
                )
                THEN 'None'
            -- not categorized: other
            ELSE NULL END AS ventilation_status
    FROM tm
    LEFT JOIN b.mimiciv_derived.ventilator_setting vs
        ON tm.stay_id = vs.stay_id
            AND tm.charttime = vs.charttime
    LEFT JOIN b.mimiciv_derived.oxygen_delivery od
        ON tm.stay_id = od.stay_id
            AND tm.charttime = od.charttime
)

 select c.stay_id,v.charttime,v.ventilation_status from r2 c left join lateral(select * from vs where vs.stay_id=c.stay_id and vs.charttime<=c.t0 and vs.charttime>c.t0-interval '14h' and vs.ventilation_status is not null order by charttime desc limit 1)v on true;
select c.stay_id,x.itemid,x.starttime,x.storetime from r2 c join b.mimiciv_icu.inputevents x using(stay_id)
    where x.itemid in(221906,222315) and x.starttime<=c.t0 and x.endtime>c.t0 and x.statusdescription is distinct from 'Rewritten';