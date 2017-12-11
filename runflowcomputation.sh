#!/bin/bash
nohup nice matlab-r2017b -nodisplay -nojvm -nosplash -nodesktop -r "try; run('runflowcomputation.m'); catch ME; disp(ME); disp(getReport(ME, 'extended', 'hyperlinks', 'off')); exit(1); end; exit(0);" > runflowcomputation.log 2>&1
