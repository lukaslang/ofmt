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
groupanalysis = true;

% Flag for region analysis.
regionanalysis = true;

% Flag for individual analysis.
individualanalysis = true;

% Create output folder.
outputFolder = fullfile(resultfolder);
mkdir(outputFolder);

% Set result folder.
resultFolder = fullfile('results', 'flow');

% Add all subfolders.
y = dir(resultFolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Combined analysis.
for k=1:length(groups)
    groupname = groups(k).name;
    groupfolder = fullfile(resultFolder, groupname);

    fprintf('Group: %s\n', groupfolder);

    % Run analysis.
    [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder);

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
end

function [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder)
%LOADDATASET Groups split sequences to datasets in one group.

% Define outlier datasets.
keySet = {'11_036'};
valueSet = {true};
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
ds = cell(0, 1);
for k=1:length(uds)
    if ~isKey(M, uds{k})
        ds{k} = uds{k};
    end
end

% Iterate through unique datasets.
for k=1:length(ds)
    fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
    [v1{k}, v2{k}, seg{k}, seg1{k}, seg2{k}, f{k}, u{k}] = loaddataset(groupfolder, ds{k});
end
end

function [v1, v2, seg, seg1, seg2, f, u] = loaddataset(groupfolder, dataset)

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
    if(~(exist(fullfile(folder, 'results-denoising.mat'), 'file') && exist(fullfile(folder, 'results-flow.mat'), 'file')))
        warning('No result or data found for sequence: %s.\n', folder);
        continue;
    end
    
    % Load data.
    F = load(fullfile(folder, 'results-denoising.mat'), 'fdelta', 'f');
    f = cat(3, f, F.fdelta);
    u = cat(3, u, F.f);
    
    V = load(fullfile(folder, 'results-flow.mat'), 'v1', 'v2');
    v1 = cat(3, v1, V.v1);
    v2 = cat(3, v2, V.v2);

    % Load segmentation of first sequence.
    if(l == 1)
        seg = im2double(imread(fullfile(folder, 'segmentation.png')));
        if(isfile(fullfile(folder, 'segmentation1.png')) && isfile(fullfile(folder, 'segmentation2.png')))
            seg1 = im2double(imread(fullfile(folder, 'segmentation1.png')));
            seg2 = im2double(imread(fullfile(folder, 'segmentation2.png')));
        end
    end
    l = l + 1;
end

% Set pixel size.
keySet = {'02_006', '02_009', '02_011', '02_013', '04_002', '04_004', '04_010', '04_015',...
          '16_005', '16_009', '17_010', '17_014', '17_025', '17_028', '17_031', '17_034', '17_040',...
          '01_011', '04_005', '04_007', '04_011', '06_005', '06_009', '06_012', '06_015', '06_022',...
          '09_014', '09_022', '09_024',...
          '10_020', '10_023', '11_031', '11_036', '05_004', '05_003', '05_005', '05_008', '05_020',...
          '10_013', '10_017', '10_027', '11_004', '11_018', '14_007', '01_015', '01_017', '01_021', '01_023',...
          '12_014', '12_018', '12_022', '12_029', '12_033', '12_037', '13_008', '13_012', '13_018', '13_022'};
valueSet = [0.303, 0.217, 0.23, 0.303, 0.188, 0.257, 0.269, 0.298,...
            0.223, 0.168, 0.24, 0.276, 0.303, 0.284, 0.303, 0.286, 0.225,...
            0.217, 0.257, 0.257, 0.269, 0.281, 0.283, 0.265, 0.265, 0.228,...
            0.214, 0.255, 0.255,...
            0.335, 0.412, 0.378, 0.227, 0.22, 0.276, 0.304, 0.28, 0.304,...
            0.265, 0.335, 0.381, 0.374, 0.289, 0.37, 0.303, 0.304, 0.248, 0.244,...
            0.39, 0.297, 0.297, 0.434, 0.303, 0.307, 0.324, 0.299, 0.299, 0.299];
% Sanity check and create map.
assert(length(keySet) == length(unique(keySet)));
pixelSize = containers.Map(keySet,valueSet);

% Set interval between frames (seconds).
interval = 0.65;

% Scale velocities according to pixel size.
v1 = v1 * pixelSize(dataset) / interval;
v2 = v2 * pixelSize(dataset) / interval;

end