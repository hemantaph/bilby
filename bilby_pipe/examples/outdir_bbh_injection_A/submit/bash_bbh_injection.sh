#!/usr/bin/env bash

# bbh_injection_data0_0_generation
# PARENTS 
# CHILDREN bbh_injection_data0_0_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_generation outdir_bbh_injection_A/bbh_injection_config_complete.ini --local --label bbh_injection_data0_0_generation --idx 0 --trigger-time 0 --injection-file outdir_bbh_injection_A/data/bbh_injection_injection_file.dat

# bbh_injection_data0_0_analysis_H1L1_dynesty
# PARENTS bbh_injection_data0_0_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_analysis outdir_bbh_injection_A/bbh_injection_config_complete.ini --local --outdir outdir_bbh_injection_A --detectors H1 --detectors L1 --label bbh_injection_data0_0_analysis_H1L1_dynesty --data-dump-file outdir_bbh_injection_A/data/bbh_injection_data0_0_generation_data_dump.pickle --sampler dynesty

