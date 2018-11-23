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

% This script runs parameter studies for the joint approach.
clear;
close all;
clc;

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Define and create folder with results.
resultfolder = fullfile('results', 'parameterstudies');
mkdir(resultfolder);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Set parameter range (alpha, beta, gamma) and generate combinations.
params = {[5e-3, 1e-2, 5e-2], [1e-2, 1e-1, 1e0], [1e-1, 5e-1, 1e0]};
[param1, param2, param3] = ndgrid(params{:});
ncombs = numel(param1);

% Run through all groups.
for k=1:length(groups)
    groupname = groups(k).name;
    % Run through all datasets.
    y = dir(fullfile(datapath, groupname));
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    for l=1:length(datasets)
        dataset = datasets(l).name;
        datafolder = fullfile(datapath, groupname, dataset, filesep);
        outputfolder = fullfile(resultfolder, removebrackets(groupname), dataset);
        
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        mkdir(outputfolder);
        
        % Compute and save results.
        for p=1:ncombs
            alpha = param1(p);
            beta = param2(p);
            gamma = param3(p);
            fprintf('Parameter combination %.2i/%.2i.\n', p, ncombs);
            tic;
            [f, uinit, u, v] = runjointmodel(datafolder, alpha, beta, gamma);
            toc;
            save(fullfile(outputfolder, sprintf('results-%.2i.mat', p)), 'f', 'uinit', 'u', 'v', 'alpha', 'beta', 'gamma');
        end
    end
end

function [f, uinit, u, v] = runjointmodel(datafolder, alpha, beta, gamma)
    % Load files and save mat file.
    folderContent = dir(fullfile(datafolder, '*.tif'));
    for k=1:numel(folderContent)
        rawImage = imread(fullfile(datafolder, folderContent(k).name));
        f(:, :, k) = im2double(rawImage);
    end
    
    temporalSmoothness = 0.5;
    mainJoint = jointModelLargeScale(f, alpha, beta, gamma, 'temporalSmoothness', temporalSmoothness, 'verbose', 1);
    mainJoint.doWarping = false;
    mainJoint.opticalFlowTerm = 'classic';
    mainJoint.numMainIt = 10;
    
    % Init and get denoised sequence.
    mainJoint.init; 
    uinit = mainJoint.u;

    % Run algorithm.
    mainJoint.run;

    % Recover solution.
    u = mainJoint.u;
    v = mainJoint.v;
end