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
% This script loads an image sequence, performs L2-TV-L2 denoising, and
% stores the output.
clear;
close all;
clc;

% Load data.
foldername = 'control/';
name = '02_009/';
folder = fullfile(datapath, foldername, name);
folderContent = dir(fullfile(folder, '*.tif'));

% Load files.
for k=1:numel(folderContent)
    rawImage = imread(fullfile(folder, folderContent(k).name));
    fdelta{k} = im2double(rawImage);
end

% Restrict frames.
frames = 1:5;
%frames = 1:length(f);

% Restrict image.
rngx = 200:370;
rngy = 150:256;

% Apply selection.
fdelta = cat(3, fdelta{frames});
fdelta = fdelta(rngx, rngy, :);
[n, m, t] = size(fdelta);

% Define denoising problem.
p = @denoise3dl2tv;

% Set regularisation parameters.
alpha = 0.005;
beta = 0.75;

% Set algorithm parameters.
tau = 1/sqrt(8);
sigma = 1/sqrt(8);
% Define termination criterion.
term = @(iter, p) iter > 1000;
% Define verbosity, logging, set plotting callback.
alg = pdhg(tau, sigma, term, 100, @logenergy);

% Set up problem.
p = p(fdelta, alpha, beta);

% Run algorithm.
stats = alg.run(p);

% Recover solution.
f = p.solution;

% Compute mean of first frame.
frame = f(:, :, 1);
meanf = mean(frame(:));

% Normalise image sequence so that all frames have the same mean.
fnorm = f(:, :, 1:end-1);
for k=2:t-1
    frame = f(:, :, k);
	fnorm(:, :, k) = fnorm(:, :, k) * meanf / mean(frame(:));
end

% Create output folder and save results.
outputFolder = fullfile('results', foldername, name);
mkdir(outputFolder);
save(fullfile(outputFolder, 'results-denoising.mat'), 'fdelta', 'f', 'fnorm', '-v7.3');
save(fullfile(outputFolder, 'results-denoising-params.mat'), 'rngx', 'rngy', 'alpha', 'beta', 'tau', 'sigma', 'term', 'stats', '-v7.3');

% Load, restrict, and output segmentation.
seg = imread(fullfile(folder, 'images', 'segmentationMap.png'));
seg = seg(rngx, rngy);
imwrite(seg, fullfile(outputFolder, 'segmentation.png'), 'png');

% Create approximate segmentation of ring.


% Crate output folder and save images.
outputFolder = fullfile(outputFolder, 'denoising');
mkdir(outputFolder);
for k=1:t
    imwrite(fdelta(:, :, k), fullfile(outputFolder, sprintf('input-%.3i.png', k)), 'png');
	imwrite(f(:, :, k), fullfile(outputFolder, sprintf('image-%.3i.png', k)), 'png');
end