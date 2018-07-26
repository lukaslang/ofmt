#!/bin/bash
env LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/local/cuda/lib64/libcudart.so /usr/local/cuda/lib64/libcusparse.so"
nohup nice matlab-r2016b -nodisplay -nojvm -nosplash -nodesktop -r "try; run('dataPipelineScript_flow.m'); catch ME; disp(ME); disp(getReport(ME, 'extended', 'hyperlinks', 'off')); exit(1); end; exit(0);" > dataPipelineScript_flow.log 2>&1
