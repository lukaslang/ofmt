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
% This script outputs an overlay of denoised sequence and the mean flow.
clear;
close all;
clc;

% Set result folder.
resultfolder = fullfile('results', 'figures');

% Create output folder.
outputFolder = fullfile(resultfolder);
mkdir(outputFolder);

% Set result folder.
resultFolder = fullfile('results', 'flow');

% Add all subfolders.
y = dir(resultFolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Select condition.
groups = groups(3);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Combined analysis.
for k=1:length(groups)
    groupname = groups(k).name;
    groupfolder = fullfile(resultFolder, groupname);

    fprintf('Group: %s\n', groupfolder);
    
    [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname, false);
    
    % Output overlay.
    secx = 1:512;
    secy = 129:256;
    outputoverlay(resultfolder, 'overlay', v1{5}, v2{5}, u{5}, seg{5}, secx, secy);
    
    % Output overlay.
    secx = 55:55+127;
    secy = 129:256;
    outputoverlay(resultfolder, 'overlay-detail-1', v1{5}, v2{5}, u{5}, seg{5}, secx, secy);
    
    % Output overlay.
    secx = 110:110+127;
    secy = 129:256;
    outputoverlay(resultfolder, 'overlay-detail-2', v1{5}, v2{5}, u{5}, seg{5}, secx, secy);
    
    % Output overlay.
    secx = 330:330+127;
    secy = 129:256;
    outputoverlay(resultfolder, 'overlay-detail-3', v1{5}, v2{5}, u{5}, seg{5}, secx, secy);

end