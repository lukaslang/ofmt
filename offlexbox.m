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
function offlexbox(datafolder, outputfolder)
%OFFLEXBOX Computes optial flow via FlexBox.
%
%   OFFLEXBOX(datafolder, outputfolder) takes a data folder name, an
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

% Set regularisation parameter.
alpha = 0.0005;
beta = 0.005;

% Create flexBox.
main = flexBox;
%main.params.showPrimals = 100;
main.params.tol = 1e-6;
main.params.verbose = 2;
main.params.tryCPP = 1;

for k=1:t-1
    % Add variables.
    v1var{k} = main.addPrimalVar([n, m]);
    v2var{k} = main.addPrimalVar([n, m]);
    
    % Create data term.
    term = L2opticalFlowTerm(1, f(:, :, k), f(:, :, k+1));
    main.addTerm(term, [v1var{k}, v2var{k}]);

    % Add regularisation term.
    main.addTerm(L1gradientIso(alpha, [n, m]), v1var{k});
    main.addTerm(L1gradientIso(alpha, [n, m]), v2var{k});
end

% Add temporal regularisation.
idOp = speye(prod([n, m]));
for k=1:t-2
    main.addTerm(L2operator(beta, 2, {-idOp, idOp}), [v1var{k}, v1var{k+1}]);
    main.addTerm(L2operator(beta, 2, {-idOp, idOp}), [v2var{k}, v2var{k+1}]);
end
    
% Run algorithm.
tic;
main.runAlgorithm;
toc;

v1 = cell(t-1, 1);
v2 = cell(t-1, 1);
for k=1:t-1
    v1{k} = main.getPrimal(v2var{k});
    v2{k} = main.getPrimal(v1var{k});
end

% Recover solution.
v1 = cat(3, v1{:});
v2 = cat(3, v2{:});

% Create colour representation.
col = zeros(n, m, 3, t-1, 'uint8');
for k=1:t-1
    col(:, :, :, k) = computeColour(v1(:, :, k), v2(:, :, k));
end

% Create output folder and save results.
mkdir(outputfolder);
save(fullfile(outputfolder, 'results-flexbox-flow.mat'), 'v1', 'v2', 'col', '-v7.3');
save(fullfile(outputfolder, 'results-flexbox-flow-params.mat'), 'alpha', 'beta', '-v7.3');

% Create colour representation and save images.
outputfolder = fullfile(outputfolder, 'flow-flexbox');
mkdir(outputfolder);
for k=1:t-1
	imwrite(col(:, :, :, k), fullfile(outputfolder, sprintf('flow-%.3i.png', k)), 'png');
end

% Save colourwheel.
imwrite(colourWheel, fullfile(outputfolder, 'colourwheel.png'));