#!/usr/bin/env bash

# gwosc_injection_data0_1180002601-0_generation
# PARENTS 
# CHILDREN gwosc_injection_data0_1180002601-0_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_generation outdir_bbh_gwosc_injection_B/gwosc_injection_config_complete.ini --label gwosc_injection_data0_1180002601-0_generation --idx 0 --trigger-time 1180002601.0

# gwosc_injection_data0_1180002601-0_analysis_H1L1_dynesty
# PARENTS gwosc_injection_data0_1180002601-0_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_analysis outdir_bbh_gwosc_injection_B/gwosc_injection_config_complete.ini --outdir outdir_bbh_gwosc_injection_B --detectors H1 --detectors L1 --label gwosc_injection_data0_1180002601-0_analysis_H1L1_dynesty --data-dump-file outdir_bbh_gwosc_injection_B/data/gwosc_injection_data0_1180002601-0_generation_data_dump.pickle --sampler dynesty

