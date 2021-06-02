#!/usr/bin/env bash

# time_domain_data0_1126259600-0_generation
# PARENTS 
# CHILDREN time_domain_data0_1126259600-0_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_generation outdir_time_domain/time_domain_config_complete.ini --time-domain-source-model=custom_bilby_pipe_function.my_module.my_function --label time_domain_data0_1126259600-0_generation --idx 0 --trigger-time 1126259600.0

# time_domain_data0_1126259600-0_analysis_H1L1_dynesty
# PARENTS time_domain_data0_1126259600-0_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_analysis outdir_time_domain/time_domain_config_complete.ini --time-domain-source-model=custom_bilby_pipe_function.my_module.my_function --outdir outdir_time_domain --detectors H1 --detectors L1 --label time_domain_data0_1126259600-0_analysis_H1L1_dynesty --data-dump-file outdir_time_domain/data/time_domain_data0_1126259600-0_generation_data_dump.pickle --sampler dynesty

