# Breast Cancer Survival: Clinical Data vs Gene Expression (METABRIC)

*Work in progress.*

Can gene expression predict breast cancer survival better than the clinical information doctors already have? This project compares a Cox model on clinical data only (model A) with a penalised Cox model on clinical data plus gene expression (model B), using the public METABRIC cohort (1,979 patients).

## Censoring

Survival data is often incomplete: for many patients, the exact time of the event is unknown. This is called censoring. There are three types:

| Type | What happened | What is known |
|---|---|---|
| Right censoring | Observation ended before the event | The event happens after time X, but not when |
| Interval censoring | The event happened between two check-ups | It happened between check-up A and check-up B |
| Left censoring | The event happened before observation started | It happened before time X, but not when |

A related problem is left truncation (delayed entry): a patient joins the study only some time after their risk started.

**In this project, only right censoring applies:**
- The event is death from breast cancer. 837 patients were alive at their last follow-up, and 497 died of other causes. For all of them, we only know that they did not die of breast cancer before their last follow-up.
- Kaplan-Meier curves, the log-rank test, the Cox model and the C-index all handle right censoring correctly.
- Interval censoring does not apply, because death dates are known exactly. It would apply to recurrence, which is only detected at check-ups.
- Left censoring and left truncation do not apply, because every patient's time starts at diagnosis.
