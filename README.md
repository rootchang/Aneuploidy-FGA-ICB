# TumorAneuploidyFGA
Test predictive power of tumor aneuploidy and fraction of copy-number alterations in prediction of patient prognosis after immunotherapy in different cancer types.

# Introduction
Identifying patients with low tumor mutational burden (TMB) that respond to cancer immunotherapy is desperately needed. Tumor aneuploidy measures chromosome-level copy number alteration, whereas fraction of genome encompassed by copy-number alterations (FGA) includes both chromosomal and focal copy-number events 1. Both tumor aneuploidy and FGA have been shown to play a role in cancer progression and predictive to cancer prognosis 1-3. Recently, Spurr et al. 4 analyzed data of an immunogenomic cohort (MSK-IMPACT) 5 of 1,660 advanced cancer patients from 10 different cancer types treated with immune checkpoint blockade (ICB) and found that in low-TMB (the bottom 80% within each cancer type) patients, lower tumor aneuploidy predicts significantly better pan-cancer overall survival following immunotherapy. This result is quite interesting. We are interested in two further questions: 1) Does tumor aneuploidy also predict significantly better survival of low-TMB patients following immunotherapy in individual cancer types? 2) Is FGA a better biomarker in predicting survival and response of low-TMB patients following immunotherapy given that FGA contains information of both chromosomal and focal copy number alteration? By comparing Kaplan-Meier survival curves of low-TMB patients with high versus low tumor aneuploidy score for all 10 individual cancer types, we show that low tumor aneuploidy only predicts significantly (p < 0.05) better survival for “cancer of unknown primary”. By using a different parameter value in calling FGA than that used in 4, we show that FGA is better than tumor aneuploidy in predicting prognosis of low-TMB patients following immunotherapy. In addition to predict pan-cancer level patient prognosis, FGA also significantly or marginally significantly predicts prognosis after immunotherapy for renal cell carcinoma, melanoma, and non-small cell lung cancer (NSCLC) in two MSK datasets.

# Usage
1. Download data (folder "tmb_mskcc_2018") for the Samstein et al. cohort at https://www.cbioportal.org/study/summary?id=tmb_mskcc_2018. 
2. Download aneuploidy scores ("Data.csv") for each sample at https://github.com/lfspurr/Aneuploidy-ICB. 
3. Download the GENIE v.7.1 release ("genie_data_cna_hg19.seg") at https://www.synapse.org/#!Synapse:syn7222066/wiki/405659. 
4. Download data ("41587_2021_1070_MOESM3_ESM.xlsx") for the Chowell et al. cohort at https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-021-01070-8/MediaObjects/41587_2021_1070_MOESM3_ESM.xlsx. 
5. Download this GitHub project and run the R Markdown pipeline.

# Support or Contact
changtiangen@gmail.com
