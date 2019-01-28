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
function createcirchistplots(resultfolder, groupname, ds, v1, v2, seg)
%CREATECIRCHISTPLOTS Creates CircHist plots and tables.

% Set plotting resolution.
resolution = '-r300';

% Set limits for histogram bars.
llim = -280000;
ulim = 1400000;

% Set epsilon for polar histograms.
epsilon = -inf;

% Polar histogram of all angles.
outputfolder = fullfile(resultfolder, 'circular-histogram', groupname);
mkdir(outputfolder);

% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon;

% Restrict angles.
theta = theta(idx);

% Plot angular histogram.
h = figure;
ch = CircHist(mod(rad2deg(theta), 360), 50);
ch.polarAxs.ThetaZeroLocation = 'right';
ch.fontSize = 20;
ch.polarAxs.LineWidth = 1.5;
ch.setRLim([llim, ulim]);
delete(ch.scaleBar);
delete(ch.polarAxs.Title);
export_fig(h, fullfile(outputfolder, sprintf('%s-circular-histogram-all.png', ds)), '-png', resolution, '-a1', '-transparent');
close(h);

end