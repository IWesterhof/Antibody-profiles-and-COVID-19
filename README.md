# Antibody-levels-and-COVID-19

## DOI preprint manuscript
doi: http://doi.org/xxx/](https://doi.org/10.1101/2024.08.29.24312767

## PROJECT OVERVIEW
VERDI (SARS-coV2 variants Evaluation in pRegnancy and paeDIatrics cohorts) project aims to generate improved evidence on the epidemiology, outcomes, prevention and treatment of variants of SARS-CoV-2 amongst children and pregnant women as a global response to the pandemic, involving cohort studies from diverse geographic and economic settings. 

## FUNDING
IW, GR, and PB were supported by the VERDI project (101045989), funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the Health and Digital Executive Agency. Neither the European Union nor the granting authority can be held responsible for them. 

## Methods and Data 

### Study design and study population
The VERDI-RECOVER household study was conducted in the Netherlands from January to March 2022 during the emergence of the Omicron BA.1/BA.2. The study protocol has been described elsewhere 11. In brief, households were eligible following a positive SARS-CoV-2 test result, either by rapid antigen detection assay or RT-PCR, in the household index case, and no positive test result among any of the other household members in the previous two weeks. Symptoms of household members were monitored on a daily basis by means of a symptom diary until 21 days after the last symptom onset. Study participants were asked to self-collect nose-throat swab (NTS) and saliva samples on days 0, 2, 4, 7, and 14, and upon onset of new COVID-19-related symptoms. Paired dried blood spots (DBS) were collected at enrollment and at least 4 weeks after enrollment or 10 days after the last symptom-onset. Additional data collection included COVID-19 vaccination status, SARS-CoV-2 infection history, and medical comorbidities, and was done using a questionnaire. 

### Diagnostic testing
NTS and saliva samples were tested for SARS-CoV-2 by RT-PCR, using the TaqPath™ COVID-19 ORF1, N and S gene targeted RT-PCR kit on a QUANTSTUDIO™5 (Thermofisher Scientific). 

To investigate the presence of antibodies, DBS were tested in a final dilution of 1:40 for immunoglobin G antibodies reactive with the target antigens using a multiplex protein microarray. Antibody binding was measured for the following human coronavirus targets: NL63 nucleocapsid protein (NP), S1 subunit, and S ectodomain; 229E NP and S1; HKU NP and S1; OC43 S ectodomain; as well as SARS-CoV-2 S1, NP, and S ectodomain. Recombinant spike proteins of S1 or the S ectodomain were expressed in HEK293 cells NP proteins were produced in Escherichia coli (for hCoV-OC43, hCoV-229E, hCoV-NL63, and hCoV-HKU1; Medix Biochemica, Espoo, Finland) or insect cells and baculovirus (SARS; Sino Biological, Beijing, China; Okba et al. 2020). More details about the laboratory process are given elsewhere (Sikkema et al., 2023) .

References: 
> Okba NMA, Müller MA, Li W, Wang C, GeurtsvanKessel CH, Corman VM et al. Severe Acute Respiratory Syndrome Coronavirus 2−Specific Antibody Responses in Coronavirus Disease Patients. Emerg Infect Dis 2020; 26: 1478–1488.

> Sikkema RS, Bruin E de, Ramakers C, Bentvelsen R, Li W, Bosch B-J et al. Reduced Seasonal Coronavirus Antibody Responses in Children Following COVID-19 Mitigation Measures, The Netherlands. Viruses 2023; 15: 212.


### Definitions 
SARS-CoV-2 infection in household members was defined as a positive RT-PCR result, with a cycle threshold (Ct) value 40 and lower, in at least one of the serial NTS or saliva samples. A household member was identified as a coprimary case when the NTS collected on the day of enrollment tested positive for SARS-CoV-2 RNA.Subsequently, a consecutive period with CT-values below 40 was defined as the duration of virus positivity. From this period, we also marked the duration of CT-values below 30 as a proxy for the infectious period. The peak in viral load was determined as the lowest detected CT-value during the infection. 

The presence of a previous SARS-CoV-2 infection prior to study enrolment was based on either a reported prior positive SARS-CoV-2 test result or on the presence of SARS NP antibodies at baseline.

Antibody status for SARS-CoV-2 and seasonal hCoV was assessed based on DBS samples. As antibody measurements based on fluorescence signals derived from the microarray on DBS samples were right censored and no protective reference levels are known, we used a non-parametric approach for analysis: for each antibody target, we grouped antibodies into ‘high’ or ‘low’ relative to the median value of the study population. Next, we calculated a cumulative score of high antibodies per coronavirus (SARS-CoV-2, hCoV-NL63, hCoV-229E, hCoV-HKU, hCoV-OC43). This cumulative antibody score per coronavirus included could vary between 0 and 3 depending on the number of recombinant viral antigens (spike ectodomain, spike S1, nucleoprotein) for which the antibody titer was above the median of the study population. For hCoV-OC43, only one antibody target was included, and no cumulative score was calculated. 

Cross-reactive antibodies are antibodies that are generated in response to one virus, such as SARS-CoV, and can also bind to antigens from another virus, such as SARS-CoV-2. Cross-protection occurs if these antibodies reduce success of subsequent infection. We use these definitions to assess how immune responses of different coronaviruses interact and contribute to the analysis of antibody levels and Omicron infection risk, symptom severity and CT-value trajectories.

SARS-CoV-2 disease severity was categorized into symptomatic disease, pauci-symptomatic, and asymptomatic episodes. Symptomatic disease was defined as: 1) onset of fever OR 2) two consecutive days with one respiratory (cough, sore throat, runny or congested nose, dyspnea) and one systemic symptom (headache, muscle ache, sweats or chills or tiredness) or with at least two respiratory symptoms. An episode was defined as pauci-symptomatic if symptoms occurred within the specified time-window but remained below the threshold for a symptomatic disease episode, and asymptomatic if no symptoms were reported. Subjects meeting the criteria for symptomatic disease additionally received a daily symptom severity score which consisted of a 5-point Likert scale per reported symptom present, except for fever, which was categorized as <38/38-39/39-40/>40 degrees Celsius. The cumulative severity score was defined as the sum of daily scores reported during the symptomatic period. We used the daily symptom data and date of positive test results to define the onset and resolution of a SARS-CoV-2 episode and to calculate the episode duration. An episode started on the day of symptom onset, which had to fall within the seven days before to seven days after the first positive test result. An episode ended on the last symptomatic day that was followed by at least two days without any symptoms. 


## Installation guide

The following R Version and packages that were used to analyze the data:

- R Version 4.4.1
> https://www.r-project.org/

- R Studio Version 2024.04.2+764
> https://rstudio.com/

- readr R package Version 2.1.5
> https://cran.r-project.org/web/packages/data.table/index.html

- tidyverse R package Version 2.0.0
> https://tidyverse.tidyverse.org

- magrittr R package Version 2.0.3
> https://magrittr.tidyverse.org

- dplyr R package Version 1.1.4
> https://dplyr.tidyverse.org

devtools R package Version 2.4.5
> https://devtools.r-lib.org

- reshape2 R package Version 1.4.4
> https://github.com/hadley/reshape

- ggplot2 R package Version 3.5.1
> https://ggplot2.tidyverse.org

- gtsummary R package Version 2.0.0
> https://github.com/ddsjoberg/gtsummary

- arsenal R package Version 3.6.3
> https://github.com/mayoverse/arsenal

- cowplot R package Version 1.1.3
> https://wilkelab.org/cowplot/

- PerformanceAnalytics R package Version 2.0.4
> https://github.com/braverock/PerformanceAnalytics

