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
function denoise(datafolder, outputfolder, alpha, beta, gpu)
%DENOISE Denoises an image sequence.
%
%   DENOISE(datafolder, outputfolder, alpha, beta, gpu) takes a data folder
%   name, an output folder name, regularisation parameters alpha, beta, and
%   a gpu device ID, and runs denoising, and outputs results.

% Load data.
folderContent = dir(fullfile(datafolder, '*.tif'));

% Load files.
for k=1:numel(folderContent)
    rawImage = imread(fullfile(datafolder, folderContent(k).name));
    fdelta{k} = im2double(rawImage);
end

% Select frames.
frames = 1:length(fdelta);

% Select image section.
rngx = 1:size(fdelta{1}, 1);
rngy = 1:size(fdelta{1}, 2);

% Apply selection.
fdelta = cat(3, fdelta{frames});
fdelta = fdelta(rngx, rngy, :);
[n, m, t] = size(fdelta);

% Define denoising problem.
p = @denoise3dl2tv;

% Set algorithm parameters.
tau = 1/sqrt(8);
sigma = 1/sqrt(8);

% Define termination criterion.
term = @(iter, p, pprev, tau, sigma) pdresidual(p, pprev, tau, sigma) < 1e-6;

% Define verbosity, logging, set plotting callback.
alg = pdhg(tau, sigma, 500, term, 500, @logenergy);

% Initialise dual variables.
y = zeros(n * m * t, 3);

% Create derivative operators.
[Dx, Dy, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);

% Check for GPU usage and transfer data.
if(gpu > 0)
    gpuDevice(gpu);
    Dx = gpuArray(Dx);
    Dy = gpuArray(Dy);
    Dt = gpuArray(Dt);
    fdelta = gpuArray(fdelta);
    y = gpuArray(y);
end

% Set up problem.
p = p(fdelta, alpha, beta, Dx, Dy, Dt, y);

% Run algorithm.
stats = alg.run(p);

% Recover solution.
if(gpu > 0)
    fgpu = p.solution;
    f = gather(fgpu);
    fdelta = gather(fdelta);
else
    f = p.solution;
end

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
mkdir(outputfolder);
save(fullfile(outputfolder, 'results-denoising.mat'), 'fdelta', 'f', 'fnorm', '-v7.3');
save(fullfile(outputfolder, 'results-denoising-params.mat'), 'rngx', 'rngy', 'alpha', 'beta', 'tau', 'sigma', 'term', 'stats', '-v7.3');

% Load, restrict, and output segmentation.
seg = imread(fullfile(datafolder, 'images', 'segmentationMap.png'));
seg = seg(rngx, rngy);
imwrite(seg, fullfile(outputfolder, 'segmentation.png'), 'png');

% Crate output folder and save images.
outputfolder = fullfile(outputfolder, 'denoising');
mkdir(outputfolder);
for k=1:t
    imwrite(fdelta(:, :, k), fullfile(outputfolder, sprintf('input-%.3i.png', k)), 'png');
	imwrite(f(:, :, k), fullfile(outputfolder, sprintf('image-%.3i.png', k)), 'png');
end