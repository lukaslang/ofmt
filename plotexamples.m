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

%%

% Figure with example of distributions all giving the same mean direction.
h = figure;

% Clearly right.
subplot(1, 4, 1);
theta{1} = 0.01 * randn(1000, 1);
plotpolarhistogram(theta{1});
%title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
%subplot(2, 4, 5);
%plotmeanangle(mangle);

% Slightly right.
subplot(1, 4, 2);
theta{2} = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
plotpolarhistogram(theta{2});
%title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
%subplot(2, 4, 6);
%plotmeanangle(mangle);

% Up and down.
subplot(1, 4, 3);
theta{3} = 0.01 * randn(1000, 1) + pi/2 - 1e-3;
theta{3} = [theta{3}; theta{3} + pi + 2e-3];
%[mangle, r] = meanangle(theta);
plotpolarhistogram(theta{3});
%title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
%subplot(2, 4, 7);
%plotmeanangle(mangle);

% Left and right.
subplot(1, 4, 4);
theta{4} = 0.01 * randn(1e6, 1);
theta{4} = [theta{4}; theta{4} + pi; 1 * randn(1e3, 1)];
%[mangle, r] = meanangle(theta);
plotpolarhistogram(theta{4});
%title(sprintf('var=%.3f', 1 - r), 'Interpreter', 'latex');
%subplot(2, 4, 8);
%plotmeanangle(mangle);

export_fig(h, fullfile(outputFolder, 'different-distributions.png'), '-png', '-r300', '-a1', '-transparent');

% All datasets combined.
h = figure;
plotpolarhistogram(theta);

export_fig(h, fullfile(outputFolder, 'different-distributions-combined.png'), '-png', '-r300', '-a1', '-transparent');


%%

% Figure with example of focus.
h = figure;

subplot(2, 4, 1);
%theta{1} = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
theta{1} = [2 * pi * rand(1e6, 1)];
%[mangle1, r1] = meanangle(theta1);
plotpolarhistogram(theta{1});
%title(sprintf('var=%.3f', 1 - r1), 'Interpreter', 'latex');

subplot(2, 4, 2);
theta{2}= [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
%[mangle2, r2] = meanangle(theta2);
plotpolarhistogram(theta{2});
%title(sprintf('var=%.3f', 1 - r2), 'Interpreter', 'latex');

subplot(2, 4, 3);
theta{3} = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
%[mangle3, r3] = meanangle(theta3);
plotpolarhistogram(theta{3});
%title(sprintf('var=%.3f', 1 - r3), 'Interpreter', 'latex');

subplot(2, 4, 4);
theta{4} = [2 * pi * rand(1e6, 1); 1 * randn(1e4, 1)];
%[mangle4, r4] = meanangle(theta4);
plotpolarhistogram(theta{4});
%title(sprintf('var=%.3f', 1 - r4), 'Interpreter', 'latex');

%subplot(2, 5, 5);
%plotmeanangle(mangle1);
%hold on;
%plotmeanangle(mangle2);
%plotmeanangle(mangle3);
%plotmeanangle(mangle4);

subplot(2, 4, 5);
theta{5}= 0.2 * randn(1000, 1) + pi/8;
%[mangle5, r5] = meanangle(theta5);
plotpolarhistogram(theta{5});
%title(sprintf('var=%.3f', 1 - r5), 'Interpreter', 'latex');

subplot(2, 4, 6);
theta{6} = 0.2 * randn(1000, 1) - pi/8;
%[mangle6, r6] = meanangle(theta6);
plotpolarhistogram(theta{6});
%title(sprintf('var=%.3f', 1 - r6), 'Interpreter', 'latex');

subplot(2, 4, 7);
theta{7} = 0.2 * randn(1000, 1) + pi/16;
%[mangle7, r7] = meanangle(theta7);
plotpolarhistogram(theta{7});
%title(sprintf('var=%.3f', 1 - r7), 'Interpreter', 'latex');

subplot(2, 4, 8);
theta{8} = 0.2 * randn(1000, 1) - pi/16;
%[mangle8, r8] = meanangle(theta8);
plotpolarhistogram(theta{8});
%title(sprintf('var=%.3f', 1 - r8), 'Interpreter', 'latex');

export_fig(h, fullfile(outputFolder, 'different-focus.png'), '-png', '-r300', '-a1', '-transparent');

h = figure;
subplot(2, 1, 1);
plotpolarhistogram(theta(1:4));

subplot(2, 1, 2);
plotpolarhistogram(theta(5:8));

%subplot(2, 5, 10);
%plotmeanangle(mangle5);
%hold on;
%plotmeanangle(mangle6);
%plotmeanangle(mangle7);
%plotmeanangle(mangle8);

export_fig(h, fullfile(outputFolder, 'different-focus-combined.png'), '-png', '-r300', '-a1', '-transparent');

%%


h = figure;
theta = rand(1000, 1) * 2 * pi;
theta = mod(rad2deg(theta), 360);
ch = CircHist(theta, [0, 90, 270, 360]);
ch.polarAxs.ThetaZeroLocation = 'right';

%%

theta = 0.01 * randn(1000, 1);
h = figure;
%set(gca, 'ColorOrder', 3);
%set(gca,'ColorOrderIndex', 3);
col = get(gca, 'ColorOrder');
polarhistogram(theta, 50, 'FaceColor', col(1, :), 'Normalization', 'probability', 'FaceAlpha', 0.3);

%%

h = figure;
[mangle, r] = meanangle([pi/4]);

polarplot([0, mangle], [0, 1], '-');
hold on;
%title('Mean direction for each dataset.', 'Interpreter', 'latex');
adjustfigure();
set(gca, 'RTickLabel', []);


%%

function plotpolarhistogram(theta)
    %polarhistogram(theta, 50, 'Normalization', 'probability', 'FaceAlpha', 0.3);
    if(iscell(theta))
        theta = cellfun(@(x) mod(rad2deg(x), 360), theta, 'UniformOutput', false);
    else
        theta = mod(rad2deg(theta), 360);
    end
    ch = CircHist(theta, 50);
    ch.polarAxs.ThetaZeroLocation = 'right';
end

function plotmeanangle(mangle)
    polarplot([0, mangle], [0, 1], '-');
end
