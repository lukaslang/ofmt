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
function creategroupplots(resultfolder, groupname, ds, v1, v2, seg)
%CREATEGROUPPLOTS Creates plots and figures for each group.

% Set plotting resolution.
resolution = '-r300';

% Set upper limit for polarhistogram.
upperlim = 0.1;

%% Direction of velocities.

% Polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-all');
mkdir(outputFolder);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Convert to polar coordinates.
        [theta, rho] = cart2pol(v1{l}, -v2{l});

        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
        hold on;
    end
end
rlim([0, upperlim]);
%title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

% Polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-time-averaged');
mkdir(outputFolder);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Compute mean over time.
        meanv1 = mean(v1{l}, 3);
        meanv2 = mean(v2{l}, 3);

        % Convert to polar coordinates.
        [theta, ~] = cart2pol(meanv1, -meanv2);

        % Find segmentation.
        idx = seg{l} > 0;

        polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
        hold on;
    end
end
rlim([0, upperlim]);
%title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-time-averaged.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

%% Left-right histograms.

% Left-right polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-all');
mkdir(outputFolder);
perc = zeros(length(ds), 2);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Convert to polar coordinates.
        [theta, rho] = cart2pol(v1{l}, -v2{l});

        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        hg = polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
        perc(l, :) = hg.Values;
        hold on;
    end
end
rlim([0, 1])
%title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
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
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
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
end
rlim([0, 1])
%title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

% Create table with left-right polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-time-averaged-table');
mkdir(outputFolder);
T = table(perc(:, 1), perc(:, 2), 'VariableNames', {'Left', 'Right'}, 'RowNames', ds);
writetable(T, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged.txt', removebrackets(groupname))), 'WriteRowNames', true);

%% Group plot.

% Group polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-all');
mkdir(outputFolder);
perc = zeros(length(ds), 4);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Convert to polar coordinates.
        [theta, rho] = cart2pol(v1{l}, -v2{l});

        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        hg = polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
        perc(l, :) = hg.Values;
        hold on;
    end
end
rlim([0, 1])
%title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
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
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-direction-group-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

% Group polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-time-averaged');
mkdir(outputFolder);
perc = zeros(length(ds), 4);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
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
end
rlim([0, 1])
%title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-time-averaged.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
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
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-direction-group-time-averaged.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

%% Mean velocities.

% Average velocity plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-velocity-all');
mkdir(outputFolder);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
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
end
%title('Mean velocity for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-velocity-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

%% Mean directions.

% Mean angle plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-direction-all');
mkdir(outputFolder);
h = figure(1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Convert to polar coordinates.
        [~, rho] = cart2pol(v1{l}, -v2{l});

        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        % Convert to polar coordinates.
        [theta, ~] = cart2pol(v1{l}(idx), -v2{l}(idx));

        % Compute mean angle.
        [mangle, r] = meanangle(theta);

        polarplot([0, mangle], [0, 1], '-');
        hold on;
    end
end
%title('Mean direction for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-direction-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

%% Speed.

% Boxplot with magnitudes of all velocities.
outputFolder = fullfile(resultfolder, 'boxplot-speed-all');
mkdir(outputFolder);
rhoall = [];
for l=1:length(ds)            
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        % Convert to polar coordinates.
        [~, tmprho] = cart2pol(v1{l}(idx), -v2{l}(idx));
        rhoall = [rhoall; tmprho];
    end
end
h = figure(1);
boxplot(rhoall, 'Labels', {groupname});
ylabel('$\mu$m/second', 'Interpreter', 'latex');
%title('Speed of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-all.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

% Create table with magnitudes.
outputFolder = fullfile(resultfolder, 'boxplot-speed-all-table');
mkdir(outputFolder);
summ = zeros(length(ds) + 1, 6);
for l=1:length(ds)            
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0;

        % Convert to polar coordinates.
        [~, rho] = cart2pol(v1{l}(idx), -v2{l}(idx));
        summ(l, :) = [mean(rho), quantile(rho, [0.25, 0.5, 0.75]), min(rho), max(rho)];
    end
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
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
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
end
h = figure(1);
boxplot(rhoall, 'Labels', {groupname});
ylabel('$\mu$m/second', 'Interpreter', 'latex');
%title('Speed of time-averaged velocity.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged.png', removebrackets(groupname))), '-png', resolution, '-a1', '-transparent');
close(h);

% Create table with magnitudes.
outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged-table');
mkdir(outputFolder);
summ = zeros(length(ds) + 1, 6);
for l=1:length(ds)      
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
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
end
% Add all velocities for group.
summ(length(ds) + 1, :) = [mean(rhoall), quantile(rhoall, [0.25, 0.5, 0.75]), min(rhoall), max(rhoall)];
T = table(summ(:, 1), summ(:, 2), summ(:, 3), summ(:, 4), summ(:, 5), summ(:, 6), 'VariableNames', {'Mean', 'Q25', 'Median', 'Q75', 'Min', 'Max'}, 'RowNames', {ds{:}, 'All'});
writetable(T, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged-table.txt', removebrackets(groupname))), 'WriteRowNames', true);
end