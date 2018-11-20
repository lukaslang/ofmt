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
function createstreamlines(resultfolder, groupname, dataset, v1, v2, seg, f, u)
%CREATESTREAMLINES Creates streamline plots for each dataset.

% Get size.
[m, n, t] = size(v1);

% Set face alpha for detail.
alpha = 1;

% Set colormap for streamlines.
cmaps = 'summer';

% Set line with of arrows and streamlines.
lineWidth = 1;

% Set seed points for streamlines.
[xs, ys] = meshgrid(1:15:n, 1:15:m);

% Compute average velocities within segmentation.
meanv1 = mean(v1, 3) .* seg;
meanv2 = mean(v2, 3) .* seg;

% Set parameters for streamline computation.
len = hypot(meanv1, meanv2);
lmax = max(len(:));
stepsize = 10 / lmax;
maxit = 100;

% Streamlines of mean velocities.
outputFolder = fullfile(resultfolder, 'streamlines-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
ax = axes;
imagesc(f(:, :, 1), 'AlphaData', alpha);
colormap(ax, 'gray');
axis image;
hold on;
colormap(cmaps);
streamlines2(cat(3, meanv1', meanv2'), [xs(:), ys(:)], stepsize, maxit, cmaps, lineWidth);
axis off;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-streamlines-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent', '-native');
close(h);

% Transpose and restrict flow to segmentation.
for p=1:t
    v1t(:, :, p) = (v1(:, :, p) .* seg)';
    v2t(:, :, p) = (v2(:, :, p) .* seg)';
end

% Set parameters for streamline computation.
len = hypot(v1t, v2t);
lmax = max(len(:));
stepsize = 100 / lmax;
maxit = t;

outputFolder = fullfile(resultfolder, 'streamlines', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
ax = axes;
imagesc(f(:, :, 1), 'AlphaData', alpha);
colormap(ax, 'gray');
axis image;
hold on;
colormap(cmaps);
streamlinesnonsteady2(cat(4, v1t, v2t), [xs(:), ys(:)], stepsize, maxit, cmaps, lineWidth);
axis off;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-streamlines.png', dataset)), '-png', '-q100', '-a1', '-transparent', '-native');
close(h);

end