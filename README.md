# landscape_dispersal

This repository contains all scripts and data to reproduce the analyses and figures for:

Walters, K.E., Barbour, K.M., Powers, J.M., Martiny, J.B.H. (2025) **Microbial dispersal into surface soil is limited on a meter scale**. The ISME Journal.

The final versions are in the folder "05_final_data_scripts_used_for_ISME_paper".

Older analyses have the following organization:

* "01_Raw_data" contains all the data and metadata used for this experiment.
* Those files are used in the scripts that are in "02_Scripts_to_process_data". These scripts largely just take the raw data and work it into forms that allow us to make graphs and perform statistics. For example, these scripts may clean out positive and negative controls, rarefy sequences, or calculate needed metrics.
* The output from those cleaning/processes scripts are stored in "03_Processed_data". This folder contains the cleaned data ready for analysis!
* Those analyses are done in "04_Scripts_to_make_figures_and_do_stats."


