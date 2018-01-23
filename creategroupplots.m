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
function creategroupplots(groupname, groupfolder, outputfolder)
%CREATEGROUPPLOTS Creates plots and figures for each group.
%
%   CREATEGROUPPLOTS(groupname, groupfolder, ouputfolder) takes the
%   group name, the folder where results are, and an ouput folder.

% Load datasets.
y = dir(groupfolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
datasets = y([y.isdir]);
for l=1:length(datasets)
    dataset = datasets(l).name;
    fprintf('Dataset: %s\n', fullfile(groupfolder, dataset));

    folder = fullfile(groupfolder, dataset);
    
    % Load data.
    %load(fullfile(folder, 'results-denoising.mat'), 'f');
    V = load(fullfile(folder, 'results-flow.mat'), 'v1', 'v2');
    v1{l} = V.v1;
    v2{l} = V.v2;
    
    % Load segmentation.
    seg{l} = im2double(imread(fullfile(folder, 'segmentation.png')));

    % Set segmentation to be one inside and zero outside.
    seg{l} = seg{l} > 0;
end

% Concatenate into arrays.
v1 = cat(4, v1{:});
v2 = cat(4, v2{:});
seg = cat(4, seg{:});

% Create output folder.
outputFolder = fullfile(outputfolder, groupname);
mkdir(outputFolder);

%% Visualise velocities within segmentation.

% Convert velocities polar coordinates.
[theta, rho] = cart2pol(v1 .* seg, -v2 .* seg);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Scatter plot for velocities within segmentation.
% h = figure(1);
% polarscatter(theta(idx), rho(idx), 10, '.');
% title('Velocities inside segmentation.', 'Interpreter', 'latex');
% set(gca, 'FontName', 'Helvetica' );
% set(gca, 'FontSize', 20);
% export_fig(h, fullfile(outputFolder, 'flow-all-scatter-inside.png'), '-png', '-q100', '-a1', '-transparent');
% close(h);

% Polar histogram for velocities within segmentation.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Histogram of velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-histogram-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Binary polar histogram for velocities within segmentation.
h = figure(1);
polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
title('Histogram of velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-all-histogram-polar-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

%% Visualise mean of velocities within segmentation. 

% Compute average velocities.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, rho] = cart2pol(meanv1 .* seg, -meanv2 .* seg);

% Find vectors where length is larger than epsilon.
idx = rho >= 1e-5;

% Scatter plot.
% h = figure(1);
% polarscatter(theta(idx), rho(idx), 10, '.');
% title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
% set(gca, 'FontName', 'Helvetica' );
% set(gca, 'FontSize', 20);
% export_fig(h, fullfile(outputFolder, 'flow-mean-scatter-inside.png'), '-png', '-q100', '-a1', '-transparent');
% close(h);

% Histogram plot.
h = figure(1);
polarhistogram(theta(idx), 50, 'Normalization', 'probability');
title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-histogram-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

% Binary polar histogram for velocities within segmentation.
h = figure(1);
polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
title('Histogram of velocities inside segmentation.', 'Interpreter', 'latex');
set(gca, 'FontName', 'Helvetica' );
set(gca, 'FontSize', 20);
export_fig(h, fullfile(outputFolder, 'flow-mean-histogram-polar-inside.png'), '-png', '-q100', '-a1', '-transparent');
close(h);

end