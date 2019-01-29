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
function [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname, issplit)
%LOADDATASETS Loads results that were output by runflowcomputation.m script.
%
%   [ds, v1, v2, seg, seg1, seg2, f, u] = LOADDATASETS(groupfolder, groupname, issplit)
%   takes a folder, a name, and a flag and loads all datasets within this
%   folder. If issplit is set to true then all results of the form ABC_01,
%   ABC_02, etc. are combined into one dataset.

% Load datasets.
y = dir(groupfolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
datasets = y([y.isdir]);

if(issplit)
    % Find unique datasets.
    for k=1:length(datasets)
        % Remove suffix '_k'.
        ds{k} = datasets(k).name(1:end-2);
    end

    % Get unique datasets.
    ds = unique(ds);

    % Iterate through unique datasets.
    for k=1:length(ds)
        fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
        [v1{k}, v2{k}, seg{k}, seg1{k}, seg2{k}, f{k}, u{k}] = loaddataset_split(groupfolder, groupname, ds{k});
    end
else
    % Iterate through datasets.
    ds = {datasets.name};
    for k=1:length(ds)
        fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
        [v1{k}, v2{k}, seg{k}, seg1{k}, seg2{k}, f{k}, u{k}] = loaddataset(groupfolder, groupname, ds{k});
    end
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

folder = fullfile(groupfolder, dataset);

% Check if result file exists.
if(~(exist(fullfile(folder, 'results-denoising.mat'), 'file') && exist(fullfile(folder, 'results-flow.mat'), 'file')))
    warning('No result or data found for sequence: %s.\n', folder);
    return;
end

% Load data.
F = load(fullfile(folder, 'results-denoising.mat'), 'fdelta', 'f');
f = F.fdelta;
u = F.f;

V = load(fullfile(folder, 'results-flow.mat'), 'v1', 'v2');
v1 = V.v1;
v2 = V.v2;

% Load segmentation of first sequence.
seg = im2double(imread(fullfile(folder, 'segmentation.png')));
if(isfile(fullfile(folder, 'segmentation1.png')) && isfile(fullfile(folder, 'segmentation2.png')))
    seg1 = im2double(imread(fullfile(folder, 'segmentation1.png')));
    seg2 = im2double(imread(fullfile(folder, 'segmentation2.png')));
end

% Scale velocities.
v1 = scale(groupname, dataset, v1);
v2 = scale(groupname, dataset, v2);

end

function [v1, v2, seg, seg1, seg2, f, u] = loaddataset_split(groupfolder, groupname, dataset)

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

% Scale velocities.
v1 = scale(groupname, dataset, v1);
v2 = scale(groupname, dataset, v2);

end

function v = scale(groupname, dataset, v)
    % Set interval between frames (seconds).
    interval = 0.65;

    % Get pixel size.
    ps = pixelsize(groupname, dataset);

    % Scale velocities according to pixel size.
    v = v * ps / interval;
end