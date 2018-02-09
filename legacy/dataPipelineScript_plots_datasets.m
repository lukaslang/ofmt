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
        loaddataset(groupname, groupfolder, resultfolder);
    end
end

function loaddataset(groupname, groupfolder, outputfolder)
%LOADDATASET Groups split sequences to datasets in one group.

% Create output folder.
outputFolder = fullfile(outputfolder, groupname);
mkdir(outputFolder);

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
ds = unique(ds);

% Iterate through unique datasets.
for k=1:length(ds)
    % Initialise.
    f = [];
    u = [];
    v1 = [];
    v2 = [];
    seg = [];
    
    fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
    
    % Iterate through every split.
    l = 1;
    while(exist(fullfile(groupfolder, sprintf('%s_%d', ds{k}, l)), 'dir'))
        folder = fullfile(groupfolder, sprintf('%s_%d', ds{k}, l));
        fprintf('Sequence: %s\n', folder);
        
        % Check if result file exists.
        if(~(exist(fullfile(folder, 'mat', 'results.mat'), 'file') && exist(fullfile(folder, 'mat', 'rawData.mat'), 'file')))
            warning('No result or data found for sequence: %s.\n', folder);
           continue; 
        end

        % Load data.
        W = load(fullfile(folder, 'mat', 'rawData.mat'), 'imageStack');
        f = W.imageStack;
        V = load(fullfile(folder, 'mat', 'results.mat'), 'v', 'u');
        u = V.u;
        v1 = cat(3, v1, V.v(:, :, :, 1));
        v2 = cat(3, v2, V.v(:, :, :, 2));

        % Load segmentation of first sequence.
        if(l == 1)
            seg = im2double(imread(fullfile(folder, 'images', 'segmentationMap.png')));
        end
        l = l + 1;
    end
    
    % Create plots.
    createplots(outputFolder, ds{k}, v1, v2, seg, f, u)
end
end

function createplots(outputFolder, dataset, v1, v2, seg, f, u)
%CREATEPLOTS Creates plots and figures for each group.

    % Replicate.
    segt = repmat(seg, 1, 1, size(v1, 3));
    
    % Find segmentation.
    idx = segt > 0;
    
    %% Visualise velocities within segmentation.

    % Convert velocities polar coordinates.
    [theta, rho] = cart2pol(v1, -v2);

    % Find vectors where length is larger than epsilon.
    idx = idx & abs(rho) >= 1e-3;

    % Scatter plot for velocities within segmentation.
%     h = figure(1);
%     polarscatter(theta(idx), rho(idx), 10, '.');
%     title('Velocities inside segmentation.', 'Interpreter', 'latex');
%     set(gca, 'FontName', 'Helvetica' );
%     set(gca, 'FontSize', 20);
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-scatter-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);

    % Polar histogram for velocities within segmentation.
    h = figure(1);
    polarhistogram(theta(idx), 50, 'Normalization', 'probability');
    title('Histogram of angles inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-histogram-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Binary polar histogram for velocities within segmentation.
%     h = figure(1);
%     polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
%     title('Histogram of angles inside segmentation.', 'Interpreter', 'latex');
%     set(gca, 'FontName', 'Helvetica' );
%     set(gca, 'FontSize', 20);
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-histogram-polar-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);

    %% Visualise mean of velocities within segmentation.

    % Compute average velocities within segmentation.
    meanv1 = mean(v1, 3) .* seg;
    meanv2 = mean(-v2, 3) .* seg;

    % Convert to polar coordinates.
    [theta, rho] = cart2pol(meanv1, meanv2);
    
    % Find segmentation.
    idx = seg > 0;
    
    % Find vectors where length is larger than epsilon.
    idx = idx & abs(rho) >= 1e-3;

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
%     h = figure(1);
%     polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability');
%     title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
%     set(gca, 'FontName', 'Helvetica' );
%     set(gca, 'FontSize', 20);
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-polar-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);
    
    % Colour-coding of mean velocities.
    h = figure(1);
    imagesc(flowToColorV2(cat(3, meanv1, -meanv2)));
    axis image;
    title('Colour-coding of mean velocities inside segmentation.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % First frame of noisy image.
    h = figure(1);
    colormap gray;
    imagesc(f(:, :, 1));
    daspect([1, 1, 1]);
    axis off;
    title('First frame of noisy input.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % First frame of reconstructed image.
    h = figure(1);
    colormap gray;
    imagesc(u(:, :, 1));
    daspect([1, 1, 1]);
    axis off;
    title('First frame of reconstruction.', 'Interpreter', 'latex');
    set(gca, 'FontName', 'Helvetica' );
    set(gca, 'FontSize', 20);
    export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame-reconstructed.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
end