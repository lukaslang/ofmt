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
% This script creates figures for the paper.
clear;
close all;
clc;

% Set result folder.
resultfolder = fullfile('results', 'figures');

% Flag to skip group analysis.
groupanalysis = false;

% Flag for region analysis.
regionanalysis = false;

% Flag for individual analysis.
individualanalysis = false;

% Flag to output noisy and reconstructed images, and the flow.
sequenceanalysis = false;

% Flag to output plots for thresholded images.
thresholdedanalysis = false;

% Flag for streamlines.
streamlines = false;

% Flag for control region.
controlregion = false;

% Create output folder.
outputFolder = fullfile(resultfolder);
mkdir(outputFolder);

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Combined analysis.
for k=1:length(groups)
    groupname = groups(k).name;
    groupfolder = fullfile(datapath, groupname);

    fprintf('Group: %s\n', groupfolder);

    % Run analysis.
    [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname);

    % Create plots.
    if(groupanalysis)
        creategroupplots(fullfile(resultfolder, 'region-all'), groupname, ds, v1, v2, seg);
    end
    if(regionanalysis)
        if(any(~cellfun(@isempty, seg1)))
            creategroupplots(fullfile(resultfolder, 'region-1'), groupname, ds, v1, v2, seg1);
        end
        if(any(~cellfun(@isempty, seg2)))
            creategroupplots(fullfile(resultfolder, 'region-2'), groupname, ds, v1, v2, seg2);
        end
    end
    if(individualanalysis)
        % Run through datasets.
        for l=1:length(ds)
            % Create plots.
            createplots(fullfile(resultfolder, 'invididual'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
    if(sequenceanalysis)
        for l=1:length(ds)
            outputsequence(fullfile(resultfolder, 'sequences'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
    if(thresholdedanalysis)
        for l=1:length(ds)
            createthresholdedplots(fullfile(resultfolder, 'invididual'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
    if(streamlines)
        for l=1:length(ds)
            createstreamlines(fullfile(resultfolder, 'invididual'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
    if(controlregion)
        for l=1:length(ds)
            % Select region.
            X = 1:512;
            Y = 128:256;
            segr = ~(seg{l}(X, Y) > 0);
            fr = f{l}(X, Y, :);
            ur = u{l}(X, Y, :);
            v1r = v1{l}(X, Y, :);
            v2r = v2{l}(X, Y, :);
            createplots(fullfile(resultfolder, 'controlregion'), groupname, ds{l}, v1r, v2r, segr, fr, ur);
        end
    end
end

function [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname)
%LOADDATASET Groups split sequences to datasets in one group.

% Define outlier datasets.
keySet = {'11_036', '12_037'};
valueSet = {true, true};
M = containers.Map(keySet, valueSet);

% Load datasets.
y = dir(groupfolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
datasets = y([y.isdir]);

% Find unique datasets.
for k=1:length(datasets)
    % Remove suffix '_k'.
    ds{k} = datasets(k).name(1:end-2);
end

% Get unique datasets.
uds = unique(ds);

% Remove outliers.
ds = uds(cellfun(@(x) ~isKey(M, x), uds, 'UniformOutput', true));

% Iterate through unique datasets.
for k=1:length(ds)
    fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
    [v1{k}, v2{k}, seg{k}, seg1{k}, seg2{k}, f{k}, u{k}] = loaddataset(groupfolder, groupname, ds{k});
end
end

function [v1, v2, seg, seg1, seg2, f, u] = loaddataset(groupfolder, groupname, dataset)

% Initialise.
f = [];
u = [];
v1 = [];
v2 = [];
seg = [];
seg1 = [];
seg2 = [];

% Iterate through every split.
l = 1;
while(exist(fullfile(groupfolder, sprintf('%s_%d', dataset, l)), 'dir'))
    folder = fullfile(groupfolder, sprintf('%s_%d', dataset, l));
    fprintf('Sequence: %s\n', folder);

    % Check if result file exists.
    if(~(exist(fullfile(folder, 'mat', 'results.mat'), 'file') && exist(fullfile(folder, 'mat', 'rawData.mat'), 'file')))
        warning('No result or data found for sequence: %s.\n', folder);
        continue;
    end

    % Load data.
    W = load(fullfile(folder, 'mat', 'rawData.mat'), 'imageStack');
    f = cat(3, f, W.imageStack);
    V = load(fullfile(folder, 'mat', 'results.mat'), 'v', 'u');
    u = cat(3, u, V.u);
    
    % Remove last frame as jointModelLargeScale returns zeros.
    v1 = cat(3, v1, V.v(:, :, 1:end-1, 1));
    v2 = cat(3, v2, V.v(:, :, 1:end-1, 2));

    % Load segmentation of first sequence.
    if(l == 1)
        seg = im2double(imread(fullfile(folder, 'images', 'segmentationMap.png')));
        if(isfile(fullfile(folder, 'images', 'segmentationMap1.png')) && isfile(fullfile(folder, 'images', 'segmentationMap2.png')))
            seg1 = im2double(imread(fullfile(folder, 'images', 'segmentationMap1.png')));
            seg2 = im2double(imread(fullfile(folder, 'images', 'segmentationMap2.png')));
        end
    end
    l = l + 1;
end

% Set interval between frames (seconds).
interval = 0.65;

% Get pixel size.
ps = pixelsize(groupname, dataset);

% Scale velocities according to pixel size.
v1 = v1 * ps / interval;
v2 = v2 * ps / interval;

end