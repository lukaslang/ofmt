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
function of(datafolder, outputfolder, alpha, beta, gpu)
%OF Computes optical flow.
%
%   OF(datafolder, outputfolder, alpha, beta, gpu) takes a data folder
%   name, an output folder name, regularisation parameters alpha, beta, and
%   a gpu device ID, and runs optical flow computation, and outputs results.

% Load data.
folderContent = dir(fullfile(datafolder, 'image-*.png'));

% Load files.
for k=1:numel(folderContent)
    rawImage = imread(fullfile(datafolder, folderContent(k).name));
    f{k} = im2double(rawImage);
end
f = cat(3, f{:});
[n, m, t] = size(f);

% Define optical flow problem.
p = @vecof3dl2tv;

% Set algorithm parameters.
tau = 1/sqrt(8);
sigma = 1/sqrt(8);

% Define termination criterion.
term = @(iter, p, pprev, tau, sigma) pdresidual(p, pprev, tau, sigma) < 1e-6;

% Define verbosity, logging, set plotting callback.
alg = pdhg(tau, sigma, 500, term, 500, @logenergy);

% Compute partial derivatives.
fvec = f(:);
[Dx, Dy, ~] = vecderiv3dc(m, n, t-1, 1, 1, 1);
fx = reshape(Dx * fvec(1:n*m*(t-1)), n, m, t-1);
fy = reshape(Dy * fvec(1:n*m*(t-1)), n, m, t-1);
[~, ~, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);
ft = reshape(Dt * fvec, n, m, t);
ft = ft(:, :, 1:end-1);

% Initialise primal and dual variables.
v = zeros(n * m * (t-1), 2);
y = zeros(n * m * (t-1), 6);

% Create derivative operators.
[Dx, Dy, Dt] = vecderiv3dfw(m, n, t - 1, 1, 1, 1);

% Check for GPU usage and transfer data.
if(gpu > 0)
    gpuDevice(gpu);
    Dx = gpuArray(Dx);
    Dy = gpuArray(Dy);
    Dt = gpuArray(Dt);
    fx = gpuArray(fx);
    fy = gpuArray(fy);
    ft = gpuArray(ft);
    v = gpuArray(v);
    y = gpuArray(y);
end

% Set up problem.
p = p(fx, fy, ft, alpha, beta, Dx, Dy, Dt, v, y);

% Run algorithm.
stats = alg.run(p);

% Recover solution.
if(gpu > 0)
    [v1gpu, v2gpu] = p.solution;
    v1 = gather(v1gpu);
    v2 = gather(v2gpu);
else
    [v1, v2] = p.solution;
end

% Create colour representation.
col = zeros(n, m, 3, t-1, 'uint8');
for k=1:t-1
    col(:, :, :, k) = computeColour(v1(:, :, k), v2(:, :, k));
end

% Create output folder and save results.
mkdir(outputfolder);
save(fullfile(outputfolder, 'results-flow.mat'), 'v1', 'v2', 'col', '-v7.3');
save(fullfile(outputfolder, 'results-flow-params.mat'), 'alpha', 'beta', 'tau', 'sigma', 'term', 'stats', '-v7.3');

% Create colour representation and save images.
outputfolder = fullfile(outputfolder, 'flow');
mkdir(outputfolder);
for k=1:t-1
	imwrite(col(:, :, :, k), fullfile(outputfolder, sprintf('flow-%.3i.png', k)), 'png');
end

% Save colourwheel.
imwrite(colourWheel, fullfile(outputfolder, 'colourwheel.png'));