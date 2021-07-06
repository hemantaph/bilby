#!/usr/bin/env bash

# eccn_data0_1126259463-9_generation
# PARENTS 
# CHILDREN eccn_data0_1126259463-9_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_generation outdir_com/eccn_config_complete.ini --label eccn_data0_1126259463-9_generation --idx 0 --trigger-time 1126259463.9

# eccn_data0_1126259463-9_analysis_H1L1_dynesty
# PARENTS eccn_data0_1126259463-9_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_analysis outdir_com/eccn_config_complete.ini --outdir outdir_com --detectors H1 --detectors L1 --label eccn_data0_1126259463-9_analysis_H1L1_dynesty --data-dump-file outdir_com/data/eccn_data0_1126259463-9_generation_data_dump.pickle --sampler dynesty

