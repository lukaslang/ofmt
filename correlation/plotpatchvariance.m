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

% Add all subfolders.
y = dir(resultFolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Select groups.
%groups = groups([1, 2, 7]);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Combined analysis.
for k=1:length(groups)
    groupname = groups(k).name;
    groupfolder = fullfile(resultFolder, groupname);

    fprintf('Group: %s\n', groupfolder);
    
    [ds, v1, v2, seg, seg1, seg2, f, u] = loaddatasets(groupfolder, groupname, false);
    numds(k) = length(ds);

    % Compute mean for each dataset.
    meanvel = @(x) mean(x, 3);
    meanv1 = cellfun(meanvel, v1, 'UniformOutput', false);
    meanv2 = cellfun(meanvel, v2, 'UniformOutput', false);

    % Define data (either frames or mean).
    x1 = meanv1;
    x2 = meanv2;
    segx = seg;

    % Create folders.
    folder = fullfile(outputFolder, 'patch-r', removebrackets(groupname));
    mkdir(folder);
    
    % Run through all datasets.
    for l=1:length(ds)
        % Pick dataset.
        f1 = x1{l};
        f2 = x2{l};
        fseg = segx{l};

        % Get size of data.
        [m, n] = size(f1);

        % Run through patch sizes.
        %b = 4:10:24;
        b = 1:5;
        for j=1:length(b)
            % Define patch size.
            s = 2*b(j) + 1;
            
            % Extract all patches of size [s, s].
            id = @(x) {x};
            p1 = nlfilter(f1, [s, s], id);
            p2 = nlfilter(f2, [s, s], id);
            p1 = p1(:);
            p2 = p2(:);

            % Get indices of pixels inside segmentation.
            idx = fseg > 0;

            % Compute length of mean resultant vector for each patch.
            img = patchr(p1, p2, m, n);
            r{l, j, k} = img(idx);

            img = img .* double(fseg > 0);
            [img, ~] = gray2ind(img);
            imwrite(img, parula, fullfile(folder, sprintf('%s-patch-r-patchsize-%i.png', ds{l}, s)));
        end
    end
end

% Compute mean.
for k=1:length(groups)
    avgr(k, :) = mean(cell2mat(r(1:numds(k), :, k)), 1);
    stdr(k, :) = std(cell2mat(r(1:numds(k), :, k)), 1); % / sqrt(numds(k));
end

% Pixel sizes of the considered patches.
ps = 2*b + 1;

% Set plotting resolution.
resolution = '-r300';

% Create figure.
h = figure(1);
cla;
hold on;
for k=1:length(groups)
    %errorbar(avgr(k, :), stdr(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 20);
    plot(avgr(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2);
end
xlim([0.5, 8]);
ylim([0.5, 1]);
legend({'control', '\emph{capu}', '\emph{khc(slow)}', '\emph{grk}', '\emph{capu,khc(slow)}', '\emph{capu,khc(slow)/+}', '\emph{khc(null)}'}, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
axis square;
xlabel('Edge length of patch (in pixels).', 'Interpreter', 'latex');
ylabel('Average length of $\bar{r}$.', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
xticks([1, 2, 3]);
xticklabels({ps(1), ps(2), ps(3)});
export_fig(h, fullfile(outputFolder, 'patch-r', 'patch-r.png'), '-png', resolution, '-a1', '-transparent');

% Compute variance.
for k=1:length(groups)
    avgcircvar(k, :) = mean(1 - cell2mat(r(1:numds(k), :, k)), 1);
    stdcircvar(k, :) = std(1 - cell2mat(r(1:numds(k), :, k)), 1); % / sqrt(numds(k));
end

% Create figure.
h = figure(2);
cla;
hold on;
for k=1:length(groups)
    %errorbar(avgcircvar(k, :), stdcircvar(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 20);
    plot(avgcircvar(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2);
end
xlim([0.5, 8]);
ylim([0, 0.5]);
legend({'control', '\emph{capu}', '\emph{khc(slow)}', '\emph{grk}', '\emph{capu,khc(slow)}', '\emph{capu,khc(slow)/+}', '\emph{khc(null)}'}, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
axis square;
xlabel('Edge length of patch (in pixels).', 'Interpreter', 'latex');
ylabel('Average circular variance $S$.', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
xticks([1, 2, 3]);
xticklabels({ps(1), ps(2), ps(3)});
export_fig(h, fullfile(outputFolder, 'patch-r', 'patch-circular-variance.png'), '-png', resolution, '-a1', '-transparent');

%% Small.
% Compute mean.
for k=1:length(groups)
    avgr(k, :) = mean(cell2mat(r(1:numds(k), :, k)), 1);
    stdr(k, :) = std(cell2mat(r(1:numds(k), :, k)), 1); % / sqrt(numds(k));
end

% Pixel sizes of the considered patches.
ps = 2*b + 1;

% Set plotting resolution.
resolution = '-r300';

% Create figure.
h = figure(1);
cla;
hold on;
for k=1:length(groups)
    %errorbar(avgr(k, :), stdr(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 20);
    plot(avgr(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2);
end
xlim([0.5, 15]);
ylim([0.5, 1]);
legend({'control', '\emph{capu}', '\emph{khc(slow)}', '\emph{grk}', '\emph{capu,khc(slow)}', '\emph{capu,khc(slow)/+}', '\emph{khc(null)}'}, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
axis square;
xlabel('Edge length of patch (in pixels).', 'Interpreter', 'latex');
ylabel('Average length of $\bar{r}$.', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
xticks([1, 2, 3, 4, 5]);
xticklabels({ps(1), ps(2), ps(3), ps(4), ps(5)});
export_fig(h, fullfile(outputFolder, 'patch-r', 'patch-r.png'), '-png', resolution, '-a1', '-transparent');

% Compute variance.
for k=1:length(groups)
    avgcircvar(k, :) = mean(1 - cell2mat(r(1:numds(k), :, k)), 1);
    stdcircvar(k, :) = std(1 - cell2mat(r(1:numds(k), :, k)), 1); % / sqrt(numds(k));
end

% Create figure.
h = figure(2);
cla;
hold on;
for k=1:length(groups)
    %errorbar(avgcircvar(k, :), stdcircvar(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2, 'CapSize', 20);
    plot(avgcircvar(k, :), '-.', 'Marker', '*', 'MarkerSize', 5, 'LineWidth', 2);
end
xlim([0.5, 15]);
ylim([0, 0.5]);
legend({'control', '\emph{capu}', '\emph{khc(slow)}', '\emph{grk}', '\emph{capu,khc(slow)}', '\emph{capu,khc(slow)/+}', '\emph{khc(null)}'}, 'Interpreter', 'latex', 'FontSize', 20, 'FontName', 'Arial');
axis square;
xlabel('Edge length of patch (in pixels).', 'Interpreter', 'latex');
ylabel('Average circular variance $S$.', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
xticks([1, 2, 3, 4, 5]);
xticklabels({ps(1), ps(2), ps(3), ps(4), ps(5)});
export_fig(h, fullfile(outputFolder, 'patch-r', 'patch-circular-variance.png'), '-png', resolution, '-a1', '-transparent');