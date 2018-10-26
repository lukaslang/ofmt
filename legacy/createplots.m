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
function createplots(resultfolder, groupname, dataset, v1, v2, seg, f, u)
%CREATEPLOTS Creates plots and figures for each dataset.

% First frame of noisy image.
outputFolder = fullfile(resultfolder, 'first-frame-noisy', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
colormap gray;
imagesc(f(:, :, 1));
daspect([1, 1, 1]);
axis off;
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
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame-reconstructed.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Visualise mean of velocities within segmentation.
% Compute average velocities within segmentation.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Colour-coding of mean velocities.
outputFolder = fullfile(resultfolder, 'velocity-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
imagesc(flowToColorV2(cat(3, meanv1 .* seg, meanv2 .* seg)));
daspect([1, 1, 1]);
axis off;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-velocity-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Visualise magnitude of time-averaged velocities.
outputFolder = fullfile(resultfolder, 'speed-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
imagesc(hypot(meanv1 .* seg, meanv2 .* seg));
daspect([1, 1, 1]);
axis off;
colorbar;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-speed-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Set epsilon for polar histograms.
epsilon = -inf;

% Polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-all', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));
% Find segmentation.
idx = segt > 0 & rho > epsilon;
polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 0.08])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-all.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, ~] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 0.08])
title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Left-right polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-all', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon;

polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-all.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Left-right polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, ~] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Group polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-all', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon;

polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-all.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Group polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-time-averaged');
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[theta, ~] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

% Plot data.
polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Average velocity plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-velocity', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0;

% Compute average velocities within segmentation.
meanv1 = mean(v1(idx));
meanv2 = mean(v2(idx));      

% Convert to polar coordinates.
[theta, rho] = cart2pol(meanv1, -meanv2);

polarplot([0, theta], [0, rho], '-');
hold on;
title('Mean velocity for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-velocity.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Mean angle plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-direction', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[~, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon;

% Convert to polar coordinates.
[theta, ~] = cart2pol(v1(idx), -v2(idx));

% Compute mean angle.
mangle = meanangle(theta);

polarplot([0, mangle], [0, 1], '-');
hold on;
title('Mean direction for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-direction.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Boxplot with magnitudes of all velocities.
outputFolder = fullfile(resultfolder, 'boxplot-speed-all', removebrackets(groupname));
mkdir(outputFolder);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0;

% Convert to polar coordinates.
[~, rhoall] = cart2pol(v1(idx), -v2(idx));
h = figure(1);
boxplot(rhoall, 'Labels', {groupname});
ylabel('$\mu$m/second', 'Interpreter', 'latex');
title('Speed of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-all.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Boxplot with magnitudes of time-averaged velocities.
outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged', removebrackets(groupname));
mkdir(outputFolder);
% Compute mean over time.
meanv1 = mean(v1, 3);
meanv2 = mean(v2, 3);

% Convert to polar coordinates.
[~, tmprho] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

% Restrict.
rhoall = tmprho(idx);
h = figure(1);
boxplot(rhoall, 'Labels', {groupname});
ylabel('$\mu$m/second', 'Interpreter', 'latex');
title('Speed of time-averaged velocity.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);
end