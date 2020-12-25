#!/usr/bin/env bash

# bbh_injection_data0_0_generation
# PARENTS 
# CHILDREN bbh_injection_data0_0_analysis_H1L1_dynesty
/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py37/bin/bilby_pipe_generation outdir_bbh_injection/bbh_injection_config_complete.ini --local --label bbh_injection_data0_0_generation --idx 0 --trigger-time 0 --injection-file outdir_bbh_injection/data/bbh_injection_injection_file.dat

# bbh_injection_data0_0_analysis_H1L1_dynesty
# PARENTS bbh_injection_data0_0_generation
# CHILDREN 
/cvmfs/oasis.opensciencegrid.org/ligo/sw/conda/envs/igwn-py37/bin/bilby_pipe_analysis outdir_bbh_injection/bbh_injection_config_complete.ini --local --outdir outdir_bbh_injection --detectors H1 --detectors L1 --label bbh_injection_data0_0_analysis_H1L1_dynesty --data-dump-file outdir_bbh_injection/data/bbh_injection_data0_0_generation_data_dump.pickle --sampler dynesty

