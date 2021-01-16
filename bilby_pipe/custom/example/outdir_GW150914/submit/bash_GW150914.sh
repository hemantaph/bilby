#!/usr/bin/env bash

# GW150914_data0_1126259462-4_generation
# PARENTS 
# CHILDREN GW150914_data0_1126259462-4_analysis_H1L1_dynesty GW150914_data0_1126259462-4_analysis_H1_dynesty GW150914_data0_1126259462-4_analysis_L1_dynesty
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_generation outdir_GW150914/GW150914_config_complete.ini --label GW150914_data0_1126259462-4_generation --idx 0 --trigger-time 1126259462.4

# GW150914_data0_1126259462-4_analysis_H1L1_dynesty
# PARENTS GW150914_data0_1126259462-4_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_analysis outdir_GW150914/GW150914_config_complete.ini --outdir outdir_GW150914 --detectors H1 --detectors L1 --label GW150914_data0_1126259462-4_analysis_H1L1_dynesty --data-dump-file outdir_GW150914/data/GW150914_data0_1126259462-4_generation_data_dump.pickle --sampler dynesty

# GW150914_data0_1126259462-4_analysis_H1_dynesty
# PARENTS GW150914_data0_1126259462-4_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_analysis outdir_GW150914/GW150914_config_complete.ini --outdir outdir_GW150914 --detectors H1 --label GW150914_data0_1126259462-4_analysis_H1_dynesty --data-dump-file outdir_GW150914/data/GW150914_data0_1126259462-4_generation_data_dump.pickle --sampler dynesty

# GW150914_data0_1126259462-4_analysis_L1_dynesty
# PARENTS GW150914_data0_1126259462-4_generation
# CHILDREN 
/home/hemantakumar.phurailatpam/anaconda3/envs/bilby_pipe/bin/bilby_pipe_analysis outdir_GW150914/GW150914_config_complete.ini --outdir outdir_GW150914 --detectors L1 --label GW150914_data0_1126259462-4_analysis_L1_dynesty --data-dump-file outdir_GW150914/data/GW150914_data0_1126259462-4_generation_data_dump.pickle --sampler dynesty

