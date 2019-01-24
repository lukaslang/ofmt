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

% Set epsilon for polar histograms.
epsilon = -inf;

% Combined polar histogram of all angles.
outputfolder = fullfile(resultfolder, 'circular-histogram');
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
ch.setRLim([-90000, 450000]);
delete(ch.scaleBar);
delete(ch.polarAxs.Title);
export_fig(h, fullfile(outputfolder, sprintf('%s-circular-histogram-all.png', groupname)), '-png', resolution, '-a1', '-transparent');
close(h);

% Create table with magnitudes.
% outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged-table');
% mkdir(outputFolder);
% summ = zeros(length(ds) + 1, 6);
% for l=1:length(ds)      
%     % Check if map is present for dataset.
%     if(~isempty(seg{l}))
%         % Compute mean over time.
%         meanv1 = mean(v1{l}, 3);
%         meanv2 = mean(v2{l}, 3);
% 
%         % Convert to polar coordinates.
%         [~, rho] = cart2pol(meanv1, -meanv2);
% 
%         % Find segmentation.
%         idx = seg{l} > 0;
% 
%         % Restrict.
%         summ(l, :) = [mean(rho(idx)), quantile(rho(idx), [0.25, 0.5, 0.75]), min(rho(idx)), max(rho(idx))];
%     end
% end
% % Add all velocities for group.
% summ(length(ds) + 1, :) = [mean(rhoall), quantile(rhoall, [0.25, 0.5, 0.75]), min(rhoall), max(rhoall)];
% T = table(summ(:, 1), summ(:, 2), summ(:, 3), summ(:, 4), summ(:, 5), summ(:, 6), 'VariableNames', {'Mean', 'Q25', 'Median', 'Q75', 'Min', 'Max'}, 'RowNames', {ds{:}, 'All'});
% writetable(T, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged-table.txt', removebrackets(groupname))), 'WriteRowNames', true);
end