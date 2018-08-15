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
groupanalysis = true;

% Flag for region analysis.
regionanalysis = true;

% Flag for individual analysis.
individualanalysis = false;

% Set epsilon for polar histograms.
epsilon = 1e-3;

%% Create aggregated plots for each group.

if(groupanalysis)
    % Create output folder.
    outputFolder = fullfile(resultfolder);
    mkdir(outputFolder);

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
        [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder);
        
        % Polar histogram of all angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-all');
        mkdir(outputFolder);
        h = figure(1);
        for l=1:length(ds)
            % Convert to polar coordinates.
            [theta, rho] = cart2pol(v1{l}, -v2{l});

            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0 & rho > epsilon;

            hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
            hold on;
        end
        rlim([0, 0.08])
        title('Directions of all velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Polar histogram of time-averaged angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-time-averaged');
        mkdir(outputFolder);
        h = figure(1);
        for l=1:length(ds)
            % Compute mean over time.
            meanv1 = mean(v1{l}, 3);
            meanv2 = mean(v2{l}, 3);
            
            % Convert to polar coordinates.
            [theta, ~] = cart2pol(meanv1, -meanv2);
            
            % Find segmentation.
            idx = seg{l} > 0;

            hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
            hold on;
        end
        rlim([0, 0.08])
        title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Left-right polar histogram of all angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-all');
        mkdir(outputFolder);
        perc = zeros(length(ds), 2);
        h = figure(1);
        for l=1:length(ds)
            % Convert to polar coordinates.
            [theta, rho] = cart2pol(v1{l}, -v2{l});

            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0 & rho > epsilon;

            hg = polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
            perc(l, :) = hg.Values;
            hold on;
        end
        rlim([0, 1])
        title('Directions of all velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Create table with left-right polar histogram of all angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-table');
        mkdir(outputFolder);
        T = table(perc(:, 1), perc(:, 2), 'VariableNames', {'Left', 'Right'}, 'RowNames', ds);
        writetable(T, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right.txt', removebrackets(groupname))), 'WriteRowNames', true);

        % Left-right polar histogram of time-averaged angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-time-averaged');
        mkdir(outputFolder);
        perc = zeros(length(ds), 2);
        h = figure(1);
        for l=1:length(ds)
            % Compute mean over time.
            meanv1 = mean(v1{l}, 3);
            meanv2 = mean(v2{l}, 3);
            
            % Convert to polar coordinates.
            [theta, ~] = cart2pol(meanv1, -meanv2);
            
            % Find segmentation.
            idx = seg{l} > 0;

            hg = polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
            perc(l, :) = hg.Values;
            hold on;
        end
        rlim([0, 1])
        title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Create table with left-right polar histogram of time-averaged angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-time-averaged-table');
        mkdir(outputFolder);
        T = table(perc(:, 1), perc(:, 2), 'VariableNames', {'Left', 'Right'}, 'RowNames', ds);
        writetable(T, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged.txt', removebrackets(groupname))), 'WriteRowNames', true);

        % Group polar histogram of all angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-all');
        mkdir(outputFolder);
        perc = zeros(length(ds), 4);
        h = figure(1);
        for l=1:length(ds)
            % Convert to polar coordinates.
            [theta, rho] = cart2pol(v1{l}, -v2{l});

            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0 & rho > epsilon;

            hg = polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
            perc(l, :) = hg.Values;
            hold on;
        end
        rlim([0, 1])
        title('Directions of all velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Create table with group polar histogram of all angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-all-table');
        mkdir(outputFolder);
        T = table(perc(:, 1), perc(:, 2), perc(:, 3), perc(:, 4), 'VariableNames', {'deg270_330', 'deg330_30', 'deg30_90', 'deg90_270'}, 'RowNames', ds);
        writetable(T, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-all.txt', removebrackets(groupname))), 'WriteRowNames', true);
        
        % Group boxplot.
        outputFolder = fullfile(resultfolder, 'boxplot-direction-group-all');
        mkdir(outputFolder);
        h = figure(1);
        boxplot(perc, 'Labels', {'270-330', '330-30', '30-90', '90-270'});
        xlabel('Region in degrees');
        ylabel('Frequency of directions of the velocity for each region');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-direction-group-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Group polar histogram of time-averaged angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-time-averaged');
        mkdir(outputFolder);
        perc = zeros(length(ds), 4);
        h = figure(1);
        for l=1:length(ds)
            % Compute mean over time.
            meanv1 = mean(v1{l}, 3);
            meanv2 = mean(v2{l}, 3);
            
            % Convert to polar coordinates.
            [theta, ~] = cart2pol(meanv1, -meanv2);
            
            % Find segmentation.
            idx = seg{l} > 0;

            % Plot data.
            hg = polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
            perc(l, :) = hg.Values;
            hold on;
        end
        rlim([0, 1])
        title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Create table with group polar histogram of time-averaged angles.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-time-averaged-table');
        mkdir(outputFolder);
        T = table(perc(:, 1), perc(:, 2), perc(:, 3), perc(:, 4), 'VariableNames', {'deg270_330', 'deg330_30', 'deg30_90', 'deg90_270'}, 'RowNames', ds);
        writetable(T, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-time-averaged.txt', removebrackets(groupname))), 'WriteRowNames', true);
        
        % Group boxplot.
        outputFolder = fullfile(resultfolder, 'boxplot-direction-group-time-averaged');
        mkdir(outputFolder);
        h = figure(1);
        boxplot(perc, 'Labels', {'270-330', '330-30', '30-90', '90-270'});
        xlabel('Region in degrees');
        ylabel('Frequency of direcions of the time-averaged velocity for each region');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-direction-group-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Average velocity plot.
        outputFolder = fullfile(resultfolder, 'polarplot-mean-velocity');
        mkdir(outputFolder);
        h = figure(1);
        for l=1:length(ds)
            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0;
            
            % Compute average velocities within segmentation.
            meanv1 = mean(v1{l}(idx));
            meanv2 = mean(v2{l}(idx));      

            % Convert to polar coordinates.
            [theta, rho] = cart2pol(meanv1, -meanv2);

            polarplot([0, theta], [0, rho], '-');
            hold on;
        end
        title('Mean velocity for each dataset.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-velocity.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Mean angle plot.
        outputFolder = fullfile(resultfolder, 'polarplot-mean-direction');
        mkdir(outputFolder);
        h = figure(1);
        for l=1:length(ds)
            % Convert to polar coordinates.
            [~, rho] = cart2pol(v1{l}, -v2{l});
            
            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0 & rho > epsilon;

            % Convert to polar coordinates.
            [theta, ~] = cart2pol(v1{l}(idx), -v2{l}(idx));

            % Compute mean angle.
            mangle = meanangle(theta);

            polarplot([0, mangle], [0, 1], '-');
            hold on;
        end
        title('Mean direction for each dataset.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-direction.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);

        % Boxplot with magnitudes of all velocities.
        outputFolder = fullfile(resultfolder, 'boxplot-speed-all');
        mkdir(outputFolder);
        rhoall = [];
        for l=1:length(ds)            
            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0;
            
            % Convert to polar coordinates.
            [~, tmprho] = cart2pol(v1{l}(idx), -v2{l}(idx));
            rhoall = [rhoall; tmprho];
        end
        h = figure(1);
        boxplot(rhoall, 'Labels', {groupname});
        ylabel('$\mu$m/second', 'Interpreter', 'latex');
        title('Speed of all velocities.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Create table with magnitudes.
        outputFolder = fullfile(resultfolder, 'boxplot-speed-all-table');
        mkdir(outputFolder);
        summ = zeros(length(ds) + 1, 6);
        for l=1:length(ds)            
            % Replicate segmentation.
            segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

            % Find segmentation.
            idx = segt > 0;
            
            % Convert to polar coordinates.
            [~, rho] = cart2pol(v1{l}(idx), -v2{l}(idx));
            summ(l, :) = [mean(rho), quantile(rho, [0.25, 0.5, 0.75]), min(rho), max(rho)];
        end
        % Add all velocities for group.
        summ(length(ds) + 1, :) = [mean(rhoall), quantile(rhoall, [0.25, 0.5, 0.75]), min(rhoall), max(rhoall)];
        T = table(summ(:, 1), summ(:, 2), summ(:, 3), summ(:, 4), summ(:, 5), summ(:, 6), 'VariableNames', {'Mean', 'Q25', 'Median', 'Q75', 'Min', 'Max'}, 'RowNames', {ds{:}, 'All'});
        writetable(T, fullfile(outputFolder, sprintf('%s-boxplot-speed-all-table.txt', removebrackets(groupname))), 'WriteRowNames', true);
        
        % Boxplot with magnitudes of time-averaged velocities.
        outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged');
        mkdir(outputFolder);
        rhoall = [];
        for l=1:length(ds)            
            % Compute mean over time.
            meanv1 = mean(v1{l}, 3);
            meanv2 = mean(v2{l}, 3);
            
            % Convert to polar coordinates.
            [~, tmprho] = cart2pol(meanv1, -meanv2);
            
            % Find segmentation.
            idx = seg{l} > 0;
            
            % Restrict.
            rhoall = [rhoall; tmprho(idx)];
        end
        h = figure(1);
        boxplot(rhoall, 'Labels', {groupname});
        ylabel('$\mu$m/second', 'Interpreter', 'latex');
        title('Speed of time-averaged velocity.', 'Interpreter', 'latex');
        adjustfigure();
        export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
        close(h);
        
        % Create table with magnitudes.
        outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged-table');
        mkdir(outputFolder);
        summ = zeros(length(ds) + 1, 6);
        for l=1:length(ds)            
            % Compute mean over time.
            meanv1 = mean(v1{l}, 3);
            meanv2 = mean(v2{l}, 3);
            
            % Convert to polar coordinates.
            [~, rho] = cart2pol(meanv1, -meanv2);
            
            % Find segmentation.
            idx = seg{l} > 0;
            
            % Restrict.
            summ(l, :) = [mean(rho(idx)), quantile(rho(idx), [0.25, 0.5, 0.75]), min(rho(idx)), max(rho(idx))];
        end
        % Add all velocities for group.
        summ(length(ds) + 1, :) = [mean(rhoall), quantile(rhoall, [0.25, 0.5, 0.75]), min(rhoall), max(rhoall)];
        T = table(summ(:, 1), summ(:, 2), summ(:, 3), summ(:, 4), summ(:, 5), summ(:, 6), 'VariableNames', {'Mean', 'Q25', 'Median', 'Q75', 'Min', 'Max'}, 'RowNames', {ds{:}, 'All'});
        writetable(T, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged-table.txt', removebrackets(groupname))), 'WriteRowNames', true);
    end
end

%% Create aggregated plots for each region.

if(regionanalysis)
    % Create output folder.
    outputFolder = fullfile(resultfolder);
    mkdir(outputFolder);

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
        [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder);
    
        % Polar histogram of all angles for region one.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-region-1-all');
        mkdir(outputFolder);
        % Check if any dataset has a region map.
        if(any(~cellfun(@isempty, seg1)))
            h = figure(1);
            for l=1:length(ds)
                % Check if region map is present.
                if(~isempty(seg1{l}))
                    % Convert to polar coordinates.
                    [theta, rho] = cart2pol(v1{l}, -v2{l});

                    % Replicate segmentation.
                    segt = repmat(seg1{l}, 1, 1, size(v1{l}, 3));

                    % Find segmentation.
                    idx = segt > 0 & rho > epsilon;

                    hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
                    hold on;
                end
            end
            rlim([0, 0.08])
            title('Directions of all velocities in region 1.', 'Interpreter', 'latex');
            adjustfigure();
            export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-region-1-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
            close(h);
        end
            
        % Polar histogram of time-averaged angles for region one.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-region-1-time-averaged');
        mkdir(outputFolder);
        % Check if any dataset has a region map.
        if(any(~cellfun(@isempty, seg1)))
            h = figure(1);
            for l=1:length(ds)
                % Check if region map is present.
                if(~isempty(seg1{l}))
                    % Compute mean over time.
                    meanv1 = mean(v1{l}, 3);
                    meanv2 = mean(v2{l}, 3);

                    % Convert to polar coordinates.
                    [theta, ~] = cart2pol(meanv1, -meanv2);

                    % Find segmentation.
                    idx = seg1{l} > 0;

                    hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
                    hold on;
                end
            end
            rlim([0, 0.08])
            title('Directions of time-averaged velocities in region 1.', 'Interpreter', 'latex');
            adjustfigure();
            export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-region-1-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
            close(h);
        end

        % Polar histogram of all angles for region two.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-region-2-all');
        mkdir(outputFolder);
        % Check if any dataset has a region map.
        if(any(~cellfun(@isempty, seg2)))
            h = figure(1);
            for l=1:length(ds)
                % Check if region map is present.
                if(~isempty(seg2{l}))
                    % Convert to polar coordinates.
                    [theta, rho] = cart2pol(v1{l}, -v2{l});

                    % Replicate segmentation.
                    segt = repmat(seg2{l}, 1, 1, size(v1{l}, 3));

                    % Find segmentation.
                    idx = segt > 0 & rho > epsilon;

                    hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
                    hold on;
                end
            end
            rlim([0, 0.08])
            title('Directions of all velocities in region 1.', 'Interpreter', 'latex');
            adjustfigure();
            export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-region-2-all.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
            close(h);
        end
        
        % Polar histogram of time-averaged angles for region two.
        outputFolder = fullfile(resultfolder, 'polarhistogram-direction-region-2-time-averaged');
        mkdir(outputFolder);
        % Check if any dataset has a region map.
        if(any(~cellfun(@isempty, seg2)))
            h = figure(1);
            for l=1:length(ds)
                % Check if region map is present.
                if(~isempty(seg2{l}))
                    % Compute mean over time.
                    meanv1 = mean(v1{l}, 3);
                    meanv2 = mean(v2{l}, 3);

                    % Convert to polar coordinates.
                    [theta, ~] = cart2pol(meanv1, -meanv2);

                    % Find segmentation.
                    idx = seg2{l} > 0;

                    hg = polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
                    hold on;
                end
            end
            rlim([0, 0.08])
            title('Directions of time-averaged velocities in region 2.', 'Interpreter', 'latex');
            adjustfigure();
            export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-region-2-time-averaged.png', removebrackets(groupname))), '-png', '-q100', '-a1', '-transparent');
            close(h);
        end
        
        
    end
end
%%

if(individualanalysis)
    % Add all subfolders.
    y = dir(datapath);
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    groups = y([y.isdir]);

    fprintf('Starting analysis of folder: %s\n', datapath);
    fprintf('Found %i groups.\n', length(groups));

    % Individual analysis.
    for k=1:length(groups)
        groupname = groups(k).name;
        groupfolder = fullfile(datapath, groupname);

        fprintf('Group: %s\n', groupfolder);

        % Run analysis.
        [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder);

        % Run through datasets.
        for l=1:length(ds)
            % Create plots.
            createplots(resultfolder, groupname, ds{l}, v1{l}, v2{l}, seg{l}, f{l}, u{l});
        end
    end
end

%%

function str = removebrackets(str)

str = regexprep(regexprep(str, '[', '_'), ']', '_');

end

function [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder)
%LOADDATASET Groups split sequences to datasets in one group.

% Define outlier datasets.
keySet = {'11_036'};
valueSet = {true};
M = containers.Map(keySet, valueSet);

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
uds = unique(ds);

% Remove outliers.
ds = cell(0, 1);
for k=1:length(uds)
    if ~isKey(M, uds{k})
        ds{k} = uds{k};
    end
end

% Iterate through unique datasets.
for k=1:length(ds)
    fprintf('Dataset: %s\n', fullfile(groupfolder, ds{k}));
    [v1{k}, v2{k}, seg{k}, seg1{k}, seg2{k}, f{k}, u{k}] = loaddataset(groupfolder, ds{k});
end
end

function [v1, v2, seg, seg1, seg2, f, u] = loaddataset(groupfolder, dataset)

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
        if(isfile(fullfile(folder, 'images', 'segmentationMap1.png')) && isfile(fullfile(folder, 'images', 'segmentationMap2.png')))
            seg1 = im2double(imread(fullfile(folder, 'images', 'segmentationMap1.png')));
            seg2 = im2double(imread(fullfile(folder, 'images', 'segmentationMap2.png')));
        end
    end
    l = l + 1;
end

% Set pixel size.
keySet = {'02_006', '02_009', '02_011', '02_013', '04_002', '04_004', '04_010', '04_015',...
          '16_005', '16_009', '17_010', '17_014', '17_025', '17_028', '17_031', '17_034', '17_040',...
          '01_011', '04_005', '04_007', '04_011', '06_005', '06_009', '06_012', '06_015', '06_022',...
          '09_014', '09_022', '09_024',...
          '10_020', '10_023', '11_031', '11_036', '05_004', '05_003', '05_005', '05_008', '05_020',...
          '10_013', '10_017', '10_027', '11_004', '11_018', '14_007', '01_015', '01_017', '01_021', '01_023',...
          '12_014', '12_018', '12_022', '12_029', '12_033', '12_037', '13_008', '13_012', '13_018', '13_022'};
valueSet = [0.303, 0.217, 0.23, 0.303, 0.188, 0.257, 0.269, 0.298,...
            0.223, 0.168, 0.24, 0.276, 0.303, 0.284, 0.303, 0.286, 0.225,...
            0.217, 0.257, 0.257, 0.269, 0.281, 0.283, 0.265, 0.265, 0.228,...
            0.214, 0.255, 0.255,...
            0.335, 0.412, 0.378, 0.227, 0.22, 0.276, 0.304, 0.28, 0.304,...
            0.265, 0.335, 0.381, 0.374, 0.289, 0.37, 0.303, 0.304, 0.248, 0.244,...
            0.39, 0.297, 0.297, 0.434, 0.303, 0.307, 0.324, 0.299, 0.299, 0.299];
% Sanity check and create map.
assert(length(keySet) == length(unique(keySet)));
pixelSize = containers.Map(keySet,valueSet);

% Set interval between frames (seconds).
interval = 0.65;

% Scale velocities according to pixel size.
v1 = v1 * pixelSize(dataset) / interval;
v2 = v2 * pixelSize(dataset) / interval;

end

function createplots(resultfolder, groupname, dataset, v1, v2, seg, f, u)
%CREATEPLOTS Creates plots and figures for each group.

    % First frame of noisy image.
    outputFolder = fullfile(resultfolder, 'first-frame-noisy', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    colormap gray;
    imagesc(f(:, :, 1));
    daspect([1, 1, 1]);
    axis off;
    %title('First frame of noisy input.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame-noisy.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % First frame of reconstructed image.
    outputFolder = fullfile(resultfolder, 'first-frame-reconstructed', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    colormap gray;
    imagesc(u(:, :, 1));
    daspect([1, 1, 1]);
    axis off;
    %title('First frame of reconstruction.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame-reconstructed.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Replicate.
%     segt = repmat(seg, 1, 1, size(v1, 3));
    
    % Find segmentation.
%     idx = segt > 0;
    
    %% Visualise velocities within segmentation.

    % Convert velocities polar coordinates.
%     [theta, rho] = cart2pol(v1, -v2);

    % Find vectors where length is larger than epsilon.
%     idx = idx & abs(rho) >= 1e-3;

    % Scatter plot for velocities within segmentation.
%     h = figure(1);
%     polarscatter(theta(idx), rho(idx), 10, '.');
%     title('Velocities inside segmentation.', 'Interpreter', 'latex');
%     set(gca, 'FontName', font);
%     set(gca, 'FontSize', fontsize);
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-scatter-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);

    % Polar histogram for velocities within segmentation.
%     h = figure(1);
%     polarhistogram(theta(idx), 50, 'Normalization', 'probability');
%     title('Histogram of angles inside segmentation.', 'Interpreter', 'latex');
%     set(gca, 'FontName', font);
%     set(gca, 'FontSize', fontsize);
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-all-histogram-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);

    %% Visualise mean of velocities within segmentation.

    % Compute average velocities within segmentation.
    meanv1 = mean(v1, 3);
    meanv2 = mean(v2, 3);

    % Convert to polar coordinates.
    [theta, rho] = cart2pol(meanv1, -meanv2);
    
    % Find vectors inside segmentation and where length is larger than epsilon.
    idx = seg > 0; % & abs(rho) >= 1e-3;

    % Histogram plot.
    outputFolder = fullfile(resultfolder, 'flow-mean-histogram-inside', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceColor', [1, 1, 1]./3, 'FaceAlpha', 0.3);
    % rlim([0, 0.08])
    %title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % left-right histogram plot.
    outputFolder = fullfile(resultfolder, 'flow-mean-histogram-inside-left-right', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    polarhistogram(theta(idx), 'BinEdges', [-pi/2, pi/2, 3*pi/2], 'Normalization', 'probability', 'FaceColor', [1, 1, 1]./3, 'FaceAlpha', 0.3);
    %title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-inside-left-right.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % Group histogram plot.
    outputFolder = fullfile(resultfolder, 'flow-mean-histogram-inside-group', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceColor', [1, 1, 1]./3, 'FaceAlpha', 0.3);
    %title('Histogram of angles of mean velocities inside segmentation.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-inside-group.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);

    % Colour-coding of mean velocities.
    outputFolder = fullfile(resultfolder, 'flow-mean-inside', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    imagesc(flowToColorV2(cat(3, meanv1 .* seg, meanv2 .* seg)));
    daspect([1, 1, 1]);
    axis off;
    %title('Colour-coding of mean velocities inside segmentation.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % Convert to polar coordinates.
    [thetam, rhom] = cart2pol(mean(meanv1(idx)), -mean(meanv2(idx)));
    
    % Scatter plot.
%     outputFolder = fullfile(resultfolder, 'flow-mean-scatter-inside', removebrackets(groupname));
%     mkdir(outputFolder);
%     h = figure(1);
%     polarscatter(theta(idx), rho(idx), 10, [1, 1, 1]./3, '.');
%     hold on;
%     polarplot([0, thetam], [0, rhom], 'r-');
%     %title('Mean velocities inside segmentation.', 'Interpreter', 'latex');
%     adjustfigure();
%     export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-scatter-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
%     close(h);
    
    % Histogram of magnitudes.
    outputFolder = fullfile(resultfolder, 'flow-mean-histogram-magnitude-inside', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    histogram(rho(idx));
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-histogram-magnitude-inside.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % Boxplot of magnitudes.
    outputFolder = fullfile(resultfolder, 'flow-mean-magnitude-inside-boxplot', removebrackets(groupname));
    mkdir(outputFolder);
    h = figure(1);
    boxplot(rho(idx), 'Labels', {dataset});
    ylabel('$\mu$m/second', 'Interpreter', 'latex');
    title('Mean magnitude of velocity.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-mean-magnitude-inside-boxplot.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
    
    % Boxplot of magnitudes.
    outputFolder = fullfile(resultfolder, 'flow-magnitude-inside-boxplot', removebrackets(groupname));
    mkdir(outputFolder);
    
    % Convert to polar coordinates.
    [~, rho] = cart2pol(v1(:), -v2(:));

    h = figure(1);
    boxplot(rho(idx), 'Labels', {dataset});
    ylabel('$\mu$m/second', 'Interpreter', 'latex');
    title('Magnitude of velocity.', 'Interpreter', 'latex');
    adjustfigure();
    export_fig(h, fullfile(outputFolder, sprintf('%s-flow-magnitude-inside-boxplot.png', dataset)), '-png', '-q100', '-a1', '-transparent');
    close(h);
end

function adjustfigure()
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 20);
end