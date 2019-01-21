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
% This script plots correlation between microtubules and vesicles.
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

% Load microtubule dataset.
groupname = '08_microtubules';
groupfolder = fullfile(resultFolder, groupname);
[dsmt, v1mt, v2mt, segmt, seg1mt, seg2mt, fmt, umt] = loaddatasets(groupfolder);

% Load vesicles dataset.
groupname = '09_vesicles';
groupfolder = fullfile(resultFolder, groupname);
[dsves, v1ves, v2ves, segves, seg1ves, seg2ves, fves, uves] = loaddatasets(groupfolder);

% Correlate frame by frame.
epsilon = 1e-1;

normcomp = @(x, y) cellfun(@(x, y) x ./ sqrt(x.^2 + y.^2 + epsilon.^2), x, y, 'UniformOutput', false);
v1mtcor = normcomp(v1mt, v2mt);
v2mtcor = normcomp(v2mt, v1mt);
v1vescor = normcomp(v1ves, v2ves);
v2vescor = normcomp(v2ves, v1ves);

%normconst = @(x, y) cellfun(@(x, y) 0 ./ sqrt(x.^2 + y.^2 + epsilon.^2), x, y, 'UniformOutput', false);
%v3mtcor = normconst(v1mt, v2mt);
%v3vescor = normconst(v1ves, v2ves);

ip = cellfun(@(x1, x2, y1, y2) x1 .* y1 + x2 .* y2, v1mtcor, v2mtcor, v1vescor, v2vescor, 'UniformOutput', false);
corr = cellfun(@(x) acos(x), ip, 'UniformOutput', false);
corr = 1 - corr{1} / pi;
%corr = corr{1} / pi;

% Create output folder.
outputFolder = fullfile(outputFolder, 'correlation-mt-ves');
mkdir(outputFolder);

for k=1:size(corr, 3)
    imwrite(gray2ind(corr(:, :, k)), parula, fullfile(outputFolder, sprintf('correlation-%.3i.png', k)));
end