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
% This script runs denoising and optical flow for each dataset.
clear;
close all;
clc;

% Set GPU device ID (0 for none).
gpu = 0;

% Set regularisation parameters for denoising problem.
alpha1 = 0.005;
beta1 = 0.75;

% Set regularisation parameter for optical flow problem.
alpha2 = 0.0005;
beta2 = 0.005;

% Flag to recompute all results.
recompute = false;

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Define and create folder with results.
resultfolder = fullfile('results', 'flow');
mkdir(resultfolder);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Output folder set to: %s\n', resultfolder);
fprintf('Found %i groups.\n', length(groups));

% Check if segmentation map is available for each dataset.
checkSegmentationMap(groups);

% Run through all groups.
for k=1:length(groups)
    groupname = groups(k).name;
    % Run through all datasets.
    y = dir(fullfile(datapath, groupname));
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    for l=1:length(datasets)
        dataset = datasets(l).name;
        datafolder = fullfile(datapath, groupname, dataset);
        outputfolder = fullfile(resultfolder, groupname, dataset);
        
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        
        % Copy segmentation maps.
        if(recompute || ~exist(outputfolder, 'dir'))
            mkdir(outputfolder);
            copySegmentationMaps(datafolder, outputfolder);
        end
        
        % Run denoising if results don't exist.
        if(recompute || ~exist(fullfile(outputfolder, 'results-denoising.mat'), 'file'))
            denoise(datafolder, outputfolder, alpha1, beta1, gpu);
        end
        
        % Run optical flow computation if results don't exist.
        if(recompute || ~exist(fullfile(outputfolder, 'results-flow.mat'), 'file'))
            of(fullfile(outputfolder, 'denoising'), outputfolder, alpha2, beta2, gpu);
        end
        
        % Run optical flow computation if results don't exist.
        %if(recompute || ~exist(fullfile(outputfolder, 'results-flexbox-flow.mat'), 'file'))
            %offlexbox(fullfile(outputfolder, 'denoising'), outputfolder);
        %end
    end
end

function checkSegmentationMap(groups)
% CHECKSEGMENTATIONMAP Runs a quick check if for every dataset a
% segmentation exists. Fails with error.
    for k=1:length(groups)
        groupname = groups(k).name;
        % Run through all datasets.
        y = dir(fullfile(datapath, groupname));
        y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
        datasets = y([y.isdir]);
        for l=1:length(datasets)
            dataset = datasets(l).name;
            datafolder = fullfile(datapath, groupname, dataset);
            if(~exist(fullfile(datafolder, 'images', 'segmentationMap.png'), 'file'))
                error('Segmentation map missing for dataset: %s\n', datafolder);
            end
        end
    end
end

function copySegmentationMaps(datafolder, outputfolder)
% COPYSEGMENTATIONMAPS Loads and copies segmentation maps to result folder.
    % Load and output segmentations.
    seg = imread(fullfile(datafolder, 'images', 'segmentationMap.png'));
    imwrite(seg, fullfile(outputfolder, 'segmentation.png'), 'png');
    
    if(isfile(fullfile(datafolder, 'images', 'segmentationMap1.png')))
        seg = imread(fullfile(datafolder, 'images', 'segmentationMap1.png'));
        imwrite(seg, fullfile(outputfolder, 'segmentation1.png'), 'png');
    end
    
    if(isfile(fullfile(datafolder, 'images', 'segmentationMap2.png')))
        seg = imread(fullfile(datafolder, 'images', 'segmentationMap2.png'));
        imwrite(seg, fullfile(outputfolder, 'segmentation2.png'), 'png');
    end
end