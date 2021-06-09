#!/usr/bin/env bash

# frequency_domain_data0_1126259600-0_generation
# PARENTS 
# CHILDREN frequency_domain_data0_1126259600-0_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_generation outdir_frequency_domain/frequency_domain_config_complete.ini --local --label frequency_domain_data0_1126259600-0_generation --idx 0 --trigger-time 1126259600.0

# frequency_domain_data0_1126259600-0_analysis_H1L1_dynesty
# PARENTS frequency_domain_data0_1126259600-0_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_analysis outdir_frequency_domain/frequency_domain_config_complete.ini --local --outdir outdir_frequency_domain --detectors H1 --detectors L1 --label frequency_domain_data0_1126259600-0_analysis_H1L1_dynesty --data-dump-file outdir_frequency_domain/data/frequency_domain_data0_1126259600-0_generation_data_dump.pickle --sampler dynesty

