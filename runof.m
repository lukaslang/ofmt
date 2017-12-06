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
% This script loads an image sequence, computes L2-TV-L2 optical flow, and
% stores the output.
clear;
close all;
clc;

% Load data.
foldername = 'control/';
name = '02_009/';
folder = fullfile('results', foldername, name, 'denoising');
folderContent = dir(fullfile(folder, 'image-*.png'));

% Load files.
for k=1:numel(folderContent)
    rawImage = imread(fullfile(folder, folderContent(k).name));
    f{k} = im2double(rawImage);
end
f = cat(3, f{:});
[n, m, t] = size(f);

% Define optical flow problem.
p = @vecof3dl2tv;

% Set regularisation parameter.
alpha = 0.0005;
beta = 0.005;

% Set algorithm parameters.
tau = 1/sqrt(8);
sigma = 1/sqrt(8);
% Define termination criterion.
term = @(iter, p) iter > 1000;
% Define logging handle.
% Define verbosity, logging, set plotting callback.
alg = pdhg(tau, sigma, term, 100, @logenergy);

% Compute partial derivatives.
fvec = f(:);
[Dx, Dy, ~] = vecderiv3dc(m, n, t-1, 1, 1, 1);
fx = reshape(Dx * fvec(1:n*m*(t-1)), n, m, t-1);
fy = reshape(Dy * fvec(1:n*m*(t-1)), n, m, t-1);
[~, ~, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);
ft = reshape(Dt * fvec, n, m, t);
ft = ft(:, :, 1:end-1);

% Set up problem.
p = p(fx, fy, ft, alpha, beta);

% Run algorithm.
stats = alg.run(p);

% Recover solution.
[v1, v2] = p.solution;

% Create colour representation.
col = zeros(n, m, 3, t-1, 'uint8');
for k=1:t-1
    col(:, :, :, k) = computeColour(v1(:, :, k), v2(:, :, k));
end

% Create output folder and save results.
outputFolder = fullfile('results', foldername, name);
mkdir(outputFolder);
save(fullfile(outputFolder, 'results-flow.mat'), 'v1', 'v2', 'col', '-v7.3');
save(fullfile(outputFolder, 'results-flow-params.mat'), 'alpha', 'beta', 'tau', 'sigma', 'term', 'stats', '-v7.3');

% Create colour representation and save images.
outputFolder = fullfile(outputFolder, 'flow');
mkdir(outputFolder);
for k=1:t-1
	imwrite(col(:, :, :, k), fullfile(outputFolder, sprintf('flow-%.3i.png', k)), 'png');
end

% Save colourwheel.
imwrite(colourWheel, fullfile(outputFolder, 'colourwheel.png'));