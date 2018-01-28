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
% This script runs analysis and creates plots for each dataset.
clear;
close all;
clc;

% Set result folder.
resultfolder = fullfile('results', 'datasets');

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
    if(~exist(fullfile(resultfolder, groupname), 'dir'))
        createplots(groupname, groupfolder, resultfolder);
    end
end

function createplots(groupname, groupfolder, outputfolder)
%CREATEPLOTS Creates plots and figures for each group.
%
%   CREATEPLOTS(groupname, groupfolder, ouputfolder) takes the
%   group name, the folder where results are, and an ouput folder.

% Create output folder.
outputFolder = fullfile(outputfolder, groupname);
mkdir(outputFolder);

% Load datasets.
y = dir(groupfolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
datasets = y([y.isdir]);
for l=1:length(datasets)
    dataset = datasets(l).name;
    fprintf('Dataset: %s\n', fullfile(groupfolder, dataset));

    folder = fullfile(groupfolder, dataset);
    
    % Check if result file exists.
    if(~exist(fullfile(folder, 'mat', 'results.mat'), 'file'))
        warning('No result found for dataset: %s.\n', folder);
       continue; 
    end
    
    % Load data.
    %load(fullfile(folder, 'results-denoising.mat'), 'f');
    V = load(fullfile(folder, 'mat', 'results.mat'), 'v');
    v1 = V.v(:, :, :, 1);
    v2 = V.v(:, :, :, 2);
    
    % Load segmentation.
    seg = im2double(imread(fullfile(folder, 'images', 'segmentationMap.png')));

    % Set segmentation to be one inside and zero outside.
    seg = seg > 0;

    %% Visualise velocities within segmentation.

    % Convert velocities polar coordinates.
    [theta, rho] = cart2pol(v1 .* seg, -v2 .* seg);

    % Find vectors where length is larger than epsilon.
    idx = abs(rho) >= 1e-3;

    % Scatter plot for velocities within segmentation.
    h = figure(1);
    polarscatter(theta(idx), rho(idx), 10, '.');
    title('Velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-scatter-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Polar histogram for velocities within segmentation.
    h = figure(1);
    polarhistogram(theta(idx), 50, 'Normalization', 'probability');
    title('Histogram of angles inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-histogram-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Binary polar histogram for velocities within segmentation.
    h = figure(1);
    polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
    title('Histogram of angles inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-histogram-polar-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    %% Visualise mean of velocities within segmentation.

    % Compute average velocities.
    meanv1 = mean(v1, 3);
    meanv2 = mean(v2, 3);

    % Convert to polar coordinates.
    [theta, rho] = cart2pol(meanv1 .* seg, -meanv2 .* seg);

    % Find vectors where length is larger than epsilon.
    idx = abs(rho) >= 1e-3;

    % Scatter plot.
    h = figure(1);
    polarscatter(theta(idx), rho(idx), 10, '.');
    title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-scatter-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Histogram plot.
    h = figure(1);
    polarhistogram(theta(idx), 50, 'Normalization', 'probability');
    title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Binary polar histogram for velocities within segmentation.
    h = figure(1);
    polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
    title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-polar-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % Colour-coding of mean velocities.
    h = figure(1);
    imagesc(flowToColorV2(cat(3, meanv1 .* seg, meanv2 .* seg)));
    axis image;
    title('Colour-coding of mean velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
end
end