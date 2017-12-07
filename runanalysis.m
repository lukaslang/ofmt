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
function runanalysis(folder)
%RUNANALYSIS Creates plots and figures.
%
%   RUNANALYSIS(folder) takes a foldername, runs analysis, and outputs
%   figures to subfolder 'analysis'.

% Load data.
load(fullfile(folder, 'results-denoising.mat'), 'f');
load(fullfile(folder, 'results-flow.mat'), 'v1', 'v2');

% Load segmentation.
seg = im2double(imread(fullfile(folder, 'segmentation.png')));

% Set segmentation to be one inside and zero outside.
seg = seg > 0;

% Create output folder.
outputFolder = fullfile(folder, 'analysis');
mkdir(outputFolder);

%% Compute mean of velocities.

% Compute means of velocities over time.
meancol = computeColour(mean(v1, 3), mean(v2, 3));

% Plot mean velocity.
h = figure(1);
imagesc(meancol);
axis image;
title('Mean velocity.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-mean.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Compute variances of velocities over time.
varv1 = var(v1, 0, 3);
varv2 = var(v2, 0, 3);

% Plot mean velocity.
h = figure(1);
imagesc(varv1);
axis image;
colorbar;
title('Variance $v_{1}$.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-variance-v1.png'), '-png', '-q100', '-a1', '-transparent');
close(h);
h = figure(1);
imagesc(varv2);
axis image;
colorbar;
title('Variance $v_{2}$.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-variance-v2.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

%% Visualise velocities within segmentation.

% Convert velocities polar coordinates.
[theta, rho] = cart2pol(v1 .* seg, -v2 .* seg);

% Scatter plot for velocities within segmentation.
h = figure(1);
polarscatter(theta(:), rho(:), 10, '.');
title('Velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-scatter-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Polar histogram for velocities within segmentation.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Histogram of velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-histogram-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

%% Visualise velocities outside segmentation.

% Convert velocities polar coordinates.
[theta, rho] = cart2pol(v1 .* ~seg, -v2 .* ~seg);

% Scatter plot for velocities within segmentation.
h = figure(1);
polarscatter(theta(:), rho(:), 10, '.');
title('Velocities outside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-scatter-outside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Polar histogram for velocities within segmentation.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Histogram of velocities outside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-histogram-outside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);


%% Visualise mean of velocities within segmentation. 

% Compute average velocities.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, rho] = cart2pol(meanv1 .* seg, -meanv2 .* seg);

% Scatter plot.
h = figure(1);
polarscatter(theta(:), rho(:), 10, '.');
title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-scatter-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Histogram plot.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-histogram-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Plot mean velocity.
h = figure(1);
imagesc(computeColour(meanv1 .* seg, -meanv2 .* seg));
axis image;
title('Mean velocity inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

%% Visualise mean of velocities outside segmentation.

% Compute average velocities.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, rho] = cart2pol(meanv1 .* ~seg, -meanv2 .* ~seg);

% Scatter plot.
h = figure(1);
polarscatter(theta(:), rho(:), 10, '.');
title('Mean velocities outside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-scatter-outside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Histogram plot.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Mean velocities outside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-histogram-outside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Plot mean velocity.
h = figure(1);
imagesc(computeColour(meanv1 .* ~seg, -meanv2 .* ~seg));
axis image;
title('Mean velocity inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-outside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

end