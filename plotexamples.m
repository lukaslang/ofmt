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
% This script just creates a few plots of distribution of directions.
clear;
close all;
clc;

outputFolder = fullfile('results', 'examples');
mkdir(outputFolder);

% Figure with example of distributions all giving the same mean direction.
h = figure;

% Clearly right.
subplot(2, 4, 1);
theta = 0.01 * randn(1000, 1);
[mangle, r] = meanangle(theta);
plotpolarhistogram(theta);
title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
subplot(2, 4, 5);
plotmeanangle(mangle);

% Slightly right.
subplot(2, 4, 2);
theta = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
[mangle, r] = meanangle(theta);
plotpolarhistogram(theta);
title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
subplot(2, 4, 6);
plotmeanangle(mangle);

% Up and down.
subplot(2, 4, 3);
theta = 0.01 * randn(1000, 1) + pi/2 - 1e-3;
theta = [theta; theta + pi + 2e-3];
[mangle, r] = meanangle(theta);
plotpolarhistogram(theta);
title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
subplot(2, 4, 7);
plotmeanangle(mangle);

% Left and right.
subplot(2, 4, 4);
theta = 0.01 * randn(1e6, 1);
theta = [theta; theta + pi; 1 * randn(1e3, 1)];
[mangle, r] = meanangle(theta);
plotpolarhistogram(theta);
title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
subplot(2, 4, 8);
plotmeanangle(mangle);

%rlim([0, 0.08])
%title('Directions of all velocities.', 'Interpreter', 'latex');
export_fig(h, fullfile(outputFolder, 'different-distributions.png'), '-png', '-r300', '-a1', '-transparent');
%close(h);

% Figure with example of focus.
h = figure;

subplot(2, 5, 1);
theta1 = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
[mangle1, r1] = meanangle(theta1);
plotpolarhistogram(theta1);
title(sprintf('var=%.3f', 1 - r1), 'Interpreter', 'latex');

subplot(2, 5, 2);
theta2 = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
[mangle2, r2] = meanangle(theta2);
plotpolarhistogram(theta2);
title(sprintf('var=%.3f', 1 - r2), 'Interpreter', 'latex');

subplot(2, 5, 3);
theta3 = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
[mangle3, r3] = meanangle(theta3);
plotpolarhistogram(theta3);
title(sprintf('var=%.3f', 1 - r3), 'Interpreter', 'latex');

subplot(2, 5, 4);
theta4 = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
[mangle4, r4] = meanangle(theta4);
plotpolarhistogram(theta4);
title(sprintf('var=%.3f', 1 - r4), 'Interpreter', 'latex');

subplot(2, 5, 5);
plotmeanangle(mangle1);
hold on;
plotmeanangle(mangle2);
plotmeanangle(mangle3);
plotmeanangle(mangle4);

subplot(2, 5, 6);
theta5 = 0.2 * randn(1000, 1) + pi/8;
[mangle5, r5] = meanangle(theta5);
plotpolarhistogram(theta5);
title(sprintf('var=%.3f', 1 - r5), 'Interpreter', 'latex');

subplot(2, 5, 7);
theta6 = 0.2 * randn(1000, 1) - pi/8;
[mangle6, r6] = meanangle(theta6);
plotpolarhistogram(theta6);
title(sprintf('var=%.3f', 1 - r6), 'Interpreter', 'latex');

subplot(2, 5, 8);
theta7 = 0.2 * randn(1000, 1) + pi/16;
[mangle7, r7] = meanangle(theta7);
plotpolarhistogram(theta7);
title(sprintf('var=%.3f', 1 - r7), 'Interpreter', 'latex');

subplot(2, 5, 9);
theta8 = 0.2 * randn(1000, 1) - pi/16;
[mangle8, r8] = meanangle(theta8);
plotpolarhistogram(theta8);
title(sprintf('var=%.3f', 1 - r8), 'Interpreter', 'latex');

subplot(2, 5, 10);
plotmeanangle(mangle5);
hold on;
plotmeanangle(mangle6);
plotmeanangle(mangle7);
plotmeanangle(mangle8);

export_fig(h, fullfile(outputFolder, 'different-focus.png'), '-png', '-r300', '-a1', '-transparent');

function plotpolarhistogram(theta)
    polarhistogram(theta, 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
end

function plotmeanangle(mangle)
    polarplot([0, mangle], [0, 1], '-');
end
