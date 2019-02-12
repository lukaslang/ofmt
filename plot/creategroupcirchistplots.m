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
function creategroupcirchistplots(resultfolder, groupname, ds, v1, v2, seg)
%CREATEGROUPCIRCHISTPLOTS Creates CircHist plots and tables.

% Set plotting resolution.
resolution = '-r300';

% Set limits for histogram bars.
llim = -280000;
ulim = 1400000;

% Set epsilon for polar histograms.
epsilon = -inf;

% Combined polar histogram of all angles.
outputfolder = fullfile(resultfolder, 'circular-histogram-all');
mkdir(outputfolder);

thetaseg = cell(length(ds), 1);
for l=1:length(ds)
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Convert to polar coordinates.
        [theta, rho] = cart2pol(v1{l}, -v2{l});

        % Replicate segmentation.
        segt = repmat(seg{l}, 1, 1, size(v1{l}, 3));

        % Find segmentation.
        idx = segt > 0 & rho > epsilon;

        % Compute angles in degrees.
        thetaseg{l} = theta(idx);
    end
end

% Plot combined angular histogram.
theta = cellfun(@(x) mod(rad2deg(x), 360), thetaseg, 'UniformOutput', false);
h = figure;
ch = CircHist(theta, 50);
ch.polarAxs.ThetaZeroLocation = 'right';
ch.fontSize = 20;
ch.polarAxs.LineWidth = 1.5;
%ch.setRLim([llim, ulim]);
delete(ch.scaleBar);
delete(ch.polarAxs.Title);
export_fig(h, fullfile(outputfolder, sprintf('%s-circular-histogram-all.png', groupname)), '-png', resolution, '-a1', '-transparent');
close(h);

% Create table with magnitudes.
outputFolder = fullfile(resultfolder, 'circular-histogram-all-table');
mkdir(outputFolder);
summ = zeros(length(ds) + 1, 5);
for l=1:length(ds)      
    % Check if map is present for dataset.
    if(~isempty(seg{l}))
        % Compute mean direction in degrees.
        mu = circ_mean(thetaseg{l});
        mu = mod(rad2deg(mu), 360);

        % Compute 95% confidence limits in degrees.
        lim = circ_confmean(thetaseg{l}, 0.05);
        lim = mod(rad2deg(lim), 360);

        % Compute length of mean vector.
        r = circ_r(thetaseg{l});

        % Compute variance.
        S = 1 - r;

        % Rayleigh test for non-uniformity.
        pval = circ_rtest(thetaseg{l});

        % Store row vector.
        summ(l, :) = [mu, lim, r, S, pval];
    end
end

% Compute quantities for combined data.
theta = cell2mat(thetaseg);

% Compute mean direction in degrees.
mu = circ_mean(theta);
mu = mod(rad2deg(mu), 360);

% Compute 95% confidence limits in degrees.
lim = circ_confmean(theta, 0.05);
lim = mod(rad2deg(lim), 360);

% Compute length of mean vector.
r = circ_r(theta);

% Compute variance.
S = 1 - r;

% Rayleigh test for non-uniformity.
pval = circ_rtest(thetaseg{l});

% Store row vector.
summ(length(ds) + 1, :) = [mu, lim, r, S, pval];

% Create and write table.
T = table(summ(:, 1), summ(:, 2), summ(:, 3), summ(:, 4), summ(:, 5), 'VariableNames', {'mean_dir', 'conf_lim', 'r', 'var', 'p_rayleigh'}, 'RowNames', {ds{:}, 'All'});
writetable(T, fullfile(outputFolder, sprintf('%s-circular-histogram-all-table.txt', removebrackets(groupname))), 'WriteRowNames', true);
end