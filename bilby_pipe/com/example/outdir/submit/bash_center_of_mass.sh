#!/usr/bin/env bash

# center_of_mass_data0_0-0_generation
# PARENTS 
# CHILDREN center_of_mass_data0_0-0_analysis_H1L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_generation outdir/center_of_mass_config_complete.ini --label center_of_mass_data0_0-0_generation --idx 0 --trigger-time 0.0

# center_of_mass_data0_0-0_analysis_H1L1_dynesty
# PARENTS center_of_mass_data0_0-0_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby/bin/bilby_pipe_analysis outdir/center_of_mass_config_complete.ini --outdir outdir --detectors H1 --detectors L1 --label center_of_mass_data0_0-0_analysis_H1L1_dynesty --data-dump-file outdir/data/center_of_mass_data0_0-0_generation_data_dump.pickle --sampler dynesty

