# Towards tailored breast cancer screening: Assessing the added value of 313-polygenic risk score and breast density in a Dutch screening cohort

BOADICEA is a multifactorial breast cancer risk prediction model that incorporates established risk factors, including polygenic risk scores (PRS) and mammographic breast density (BD). While extensively validated in high-risk populations, its performance in screening cohorts remains unclear. This study evaluates the clinical validity of BOADICEA in a population-based screening cohort.

We validated BOADICEA (v7.3.2) in the PRISMA cohort (n=82,200 women in the Dutch screening program) using a nested case-control design, including 318 cases and 771 controls. Five-year risks were calculated based on sub-models including age, “clinical” risk factors (personal hormonal, lifestyle, and disease history, along with family history up to second-degree relatives), rare pathogenic variants (PVs), BD (Volpara-based continuous and categorical), and 313-PRS. Model performance was assessed using discrimination (C-index), calibration (observed/expected ratios, slope, calibration plots), risk reclassification (Net Reclassification Index, NRI), and decision analysis (Net Benefit, NB).

# Syntax files
| File                   | Description             |
| :----                  | :----                   |
| 05_External_validation.Rmd          | Validation of the BOADICEA model in the raw dataset.
| 05_External_validation_imp.Rmd         | Validation of the BOADICEA model in an imputed dataset.
| 05_External_validation_DCIShg.Rmd       | Validation of the BOADICEA model in the dataset including high grade DCIS as cases.
| 05_External_validation_DCISall.Rmd     | Validation of the BOADICEA model in the dataset including all DCIS as cases.
| 05_External_validation_subg_Age.Rmd    | Subgroup analysis by age.
| perf_functions.R                       | Functions to calculate model performance.

# Contact
Mary Ann E. Binuya <br/>
Netherlands Cancer Institute <br/>
[m.binuya@nki.nl](m.binuya@nki.nl)

# Authors
| Author                 | Role   | Description             |
| :----                  | :----: | :----                   |
| Mary Ann Binuya   | Author | Development and support |
| Renée Verdiesen  | Author | Development and review   |
| Jim Peters    | Author | Development and review  |
| Maartje Schreurs  | Author | Review   |
| Danielle van der Waal  | Author | Review   |
| Marjanka Schmidt  | Author | Review   |
| Mireille Broeders  | Author | Review  |
