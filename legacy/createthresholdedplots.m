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
function createthresholdedplots(resultfolder, groupname, dataset, v1, v2, seg, f, u)
%CREATETHRESHOLDEDPLOTS Creates plots and figures for each dataset.

% Set threshold for comets.
threshold = 0.3;

% Find indices of pixels exceeding threshold.
idxt = u > threshold;

% Threshold reconstructed image.
ut = u;
ut(~idxt) = 0;

% Threshold velocities.
v1t = v1;
v2t = v2;
v1t(~idxt) = 0;
v2t(~idxt) = 0;

% First thresholded frame of reconstructed image.
outputFolder = fullfile(resultfolder, 'first-frame-reconstructed-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
colormap gray;
imagesc(ut(:, :, 1));
daspect([1, 1, 1]);
axis off;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-first-frame-reconstructed-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Compute average velocities within segmentation.
meanv1 = mean(v1t, 3);
meanv2 = mean(v2t, 3);

% Colour-coding of mean velocities.
outputFolder = fullfile(resultfolder, 'velocity-time-averaged-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
imagesc(flowToColorV2(cat(3, meanv1 .* seg, meanv2 .* seg)));
daspect([1, 1, 1]);
axis off;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-velocity-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Visualise magnitude of time-averaged velocities.
outputFolder = fullfile(resultfolder, 'speed-time-averaged-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
imagesc(hypot(meanv1 .* seg, meanv2 .* seg));
daspect([1, 1, 1]);
axis off;
colorbar;
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-speed-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Set epsilon for polar histograms.
epsilon = 1e-3;

% Polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-all-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));
% Find segmentation.
idx = segt > 0 & rho > epsilon & idxt;
polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 0.08])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-all-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-time-averaged-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1t, 3);
meanv2 = mean(v2t, 3);

% Convert to polar coordinates.
[theta, ~] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

polarhistogram(theta(idx), 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 0.08])
title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Left-right polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-all-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon & idxt;

polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-all-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Left-right polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-left-right-time-averaged-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1t, 3);
meanv2 = mean(v2t, 3);

% Convert to polar coordinates.
[theta, ~] = cart2pol(meanv1, -meanv2);

% Find segmentation.
idx = seg > 0;

polarhistogram(theta(idx), 'BinEdges', [pi/2, 3*pi/2, 5*pi/2], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of time-averaged velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-left-right-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Group polar histogram of all angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-all-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[theta, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon & idxt;

polarhistogram(theta(idx), 'BinEdges', [-pi/2, -pi/6, pi/6, pi/2, 9*pi/6], 'Normalization', 'probability', 'FaceAlpha', 0.3);
hold on;
rlim([0, 1])
title('Directions of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-all-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Group polar histogram of time-averaged angles.
outputFolder = fullfile(resultfolder, 'polarhistogram-direction-group-time-averaged-thresholded');
mkdir(outputFolder);
h = figure(1);
% Compute mean over time.
meanv1 = mean(v1t, 3);
meanv2 = mean(v2t, 3);

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
export_fig(h, fullfile(outputFolder, sprintf('%s-polarhistogram-direction-group-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Average velocity plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-velocity-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0;

% Compute average velocities within segmentation.
meanv1 = mean(v1t(idx));
meanv2 = mean(v2t(idx));      

% Convert to polar coordinates.
[theta, rho] = cart2pol(meanv1, -meanv2);

polarplot([0, theta], [0, rho], '-');
hold on;
title('Mean velocity for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-velocity-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Mean angle plot.
outputFolder = fullfile(resultfolder, 'polarplot-mean-direction-thresholded', removebrackets(groupname));
mkdir(outputFolder);
h = figure(1);
% Convert to polar coordinates.
[~, rho] = cart2pol(v1, -v2);

% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & rho > epsilon & idxt;

% Convert to polar coordinates.
[theta, ~] = cart2pol(v1(idx), -v2(idx));

% Compute mean angle.
mangle = meanangle(theta);

polarplot([0, mangle], [0, 1], '-');
hold on;
title('Mean direction for each dataset.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-polarplot-mean-direction-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Boxplot with magnitudes of all velocities.
outputFolder = fullfile(resultfolder, 'boxplot-speed-all-thresholded', removebrackets(groupname));
mkdir(outputFolder);
% Replicate segmentation.
segt = repmat(seg, 1, 1, size(v1, 3));

% Find segmentation.
idx = segt > 0 & idxt;

% Convert to polar coordinates.
[~, rhoall] = cart2pol(v1(idx), -v2(idx));
h = figure(1);
boxplot(rhoall, 'Labels', {groupname});
ylabel('$\mu$m/second', 'Interpreter', 'latex');
title('Speed of all velocities.', 'Interpreter', 'latex');
adjustfigure();
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-all-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);

% Boxplot with magnitudes of time-averaged velocities.
outputFolder = fullfile(resultfolder, 'boxplot-speed-time-averaged-thresholded', removebrackets(groupname));
mkdir(outputFolder);
% Compute mean over time.
meanv1 = mean(v1t, 3);
meanv2 = mean(v2t, 3);

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
export_fig(h, fullfile(outputFolder, sprintf('%s-boxplot-speed-time-averaged-thresholded.png', dataset)), '-png', '-q100', '-a1', '-transparent');
close(h);
end