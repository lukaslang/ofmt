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

% Flag for streamlines.
streamlines = false;

% Flag for control region.
controlregion = false;

% Flag for directional statistical analysis.
dirstatsanalysis = false;

% Flag for directional statistical analysis.
groupdirstatsanalysis = false;

% Set control region.
X = 1:512;
Y = 128:256;

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

    % Load datasets for group.
    [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname);

    % Create plots.
    if(groupanalysis)
        creategroupplots(fullfile(resultfolder, 'region-all'), groupname, ds, v1, v2, seg);
    end
    if(groupdirstatsanalysis)
        creategroupcirchistplots(fullfile(resultfolder, 'region-all'), groupname, ds, v1, v2, seg);
    end
    if(dirstatsanalysis)
        for l=1:length(ds)
            createcirchistplots(fullfile(resultfolder, 'individual'), groupname, ds{l}, v1{l}, v2{l}, seg{l});
        end
    end
    if(regionanalysis)
        if(any(~cellfun(@isempty, seg1)))
            creategroupplots(fullfile(resultfolder, 'region-1'), groupname, ds, v1, v2, seg1);
            creategroupcirchistplots(fullfile(resultfolder, 'region-1'), groupname, ds, v1, v2, seg1);
        end
        if(any(~cellfun(@isempty, seg2)))
            creategroupplots(fullfile(resultfolder, 'region-2'), groupname, ds, v1, v2, seg2);
            creategroupcirchistplots(fullfile(resultfolder, 'region-2'), groupname, ds, v1, v2, seg2);
        end
    end
    if(individualanalysis)
        % Run through datasets.
        for l=1:length(ds)
            % Create plots.
            createplots(fullfile(resultfolder, 'individual'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
    if(sequenceanalysis)
        for l=1:length(ds)
            outputsequence(fullfile(resultfolder, 'sequences'), groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});

            % Same for control region.
            segr = ~(seg{l}(X, Y) > 0);
            fr = f{l}(X, Y, :);
            ur = u{l}(X, Y, :);
            v1r = v1{l}(X, Y, :);
            v2r = v2{l}(X, Y, :);
            outputsequence(fullfile(resultfolder, 'sequences-controlregion'), groupname, ds{l}, v1r, v2r, segr, fr, ur);
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
            segr = ~(seg{l}(X, Y) > 0);
            fr = f{l}(X, Y, :);
            ur = u{l}(X, Y, :);
            v1r = v1{l}(X, Y, :);
            v2r = v2{l}(X, Y, :);
            createplots(fullfile(resultfolder, 'controlregion'), groupname, ds{l}, v1r, v2r, segr, fr, ur);
        end
    end
end