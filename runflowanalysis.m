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
% This script runs analysis and creates plots.
clear;
close all;
clc;

% Flag to recompute all results.
recompute = false;

% Set folder with results.
resultfolder = fullfile('results-full', 'flow');

% Set output folder.
analysisfolder = fullfile('results-full', 'analysis');

% Add all subfolders.
y = dir(resultfolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Output folder set to: %s\n', resultfolder);
fprintf('Found %i groups.\n', length(groups));

% Run through all groups.
for k=1:length(groups)
    groupname = groups(k).name;
    % Run through all datasets.
    y = dir(fullfile(resultfolder, groupname));
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    for l=1:length(datasets)
        dataset = datasets(l).name;
        outputfolder = fullfile(resultfolder, groupname, dataset);
        
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        
        % Run analysis.
        if(recompute || ~exist(fullfile(outputfolder, 'analysis'), 'dir'))
            createplots(outputfolder);
        end
    end
end

% Combined analysis.
for k=1:length(groups)
    groupname = groups(k).name;
    groupfolder = fullfile(resultfolder, groupname);
    
    fprintf('Group: %s\n', groupfolder);
    
    % Run analysis.
    if(recompute || ~exist(fullfile(analysisfolder, groupname), 'dir'))
        creategroupplots(groupname, groupfolder, analysisfolder);
    end
end