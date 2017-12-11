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
function of(datafolder, outputfolder)
%OF Computes optical flow.
%
%   OF(datafolder, outputfolder) takes a data folder name, an
%   output folder name, runs optical flow computation, and outputs results.

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