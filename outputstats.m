% Copyright 2017 Lukas Lang
%
% This file is part of OFMT.
%
%    OFMT is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    OFMT is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with OFMT.  If not, see <http://www.gnu.org/licenses/>.
%
% After running runflowcomputation.sh this script takes the
% runflowcomputation.log file and outputs some stats.
clear;
close all;
clc;

% Set file.
filename = 'runflowcomputation.log';
file = fullfile(filename);

% Load file.
fprintf('Starting analysis of file: %s\n', file);
fileID = fopen(file);
str = textscan(fileID, '%s', 'delimiter', '\n');
fclose(fileID);

% Get info about datasets.
str = str{1};
groups = matchpattern(str, '^Dataset: (.+)/.+$');
datasets = matchpattern(str, '^Dataset: .+/(.+)$');
fprintf('Number of groups processed: %i\n', length(unique(groups)));
fprintf('Number of datasets processed: %i\n', length(datasets));

% Get elapsed time.
data = matchpattern(str, '^Elapsed time is (\d+\.\d+) seconds.$');
elapsed = cellfun(@(x) str2double(x), data, 'UniformOutput', true);
elapsedden = elapsed(1:2:end);
elapsedof = elapsed(2:2:end);

% Get number of iterations.
data = matchpattern(str, '^Terminated after (\d+) iterations.$');
iter = cellfun(@(x) str2double(x), data, 'UniformOutput', true);
iterden = iter(1:2:end);
iterof = iter(2:2:end);

% Statistical analysis.
data = iterden;
fprintf('Iterations denoising:\t mean=%g, ci=%g.\n', mean(data), confint(data));

data = iterof;
fprintf('Iterations optical flow: mean=%g, ci=%g.\n', mean(data), confint(data));

data = elapsedden / 60;
fprintf('Time denoising:\t\t mean=%g min., ci=%g min.\n', mean(data), confint(data));

data = elapsedof / 60;
fprintf('Time optical flow:\t mean=%g min., ci=%g min.\n', mean(data), confint(data));

data = (elapsedden + elapsedof) / 60;
fprintf('Time both:\t\t mean=%g min., ci=%g min.\n', mean(data), confint(data));
fprintf('Time total:\t\t %g h.\n', sum(data) / 60);

function ci = confint(data)
    % Compute confidence interval.
    n = length(data);
    ts = tinv(0.975, n - 1);
    ci = ts * std(data) / sqrt(n);
end

function data = matchpattern(str, pattern)
    data = regexp(str, pattern, 'tokens', 'once');
    idx = cellfun(@(x) ~isempty(x), data, 'UniformOutput', true);
    data = cellfun(@(x) x{:}, data(idx), 'UniformOutput', false);
end
