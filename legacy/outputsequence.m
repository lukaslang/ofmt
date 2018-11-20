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
function outputsequence(resultfolder, groupname, dataset, v1, v2, seg, f, u)
%OUTPUTSEQUENCE Outputs original images, denoised images, and the flow.

% Noisy images.
outputFolder = fullfile(resultfolder, removebrackets(groupname), removebrackets(dataset), 'noisy');
mkdir(outputFolder);
for k=1:size(f, 3)
    imwrite(f(:, :, k), fullfile(outputFolder, sprintf('%s-noisy-%.3i.png', dataset, k)));
end

% Denoised images.
outputFolder = fullfile(resultfolder, removebrackets(groupname), removebrackets(dataset), 'denoised');
mkdir(outputFolder);
for k=1:size(f, 3)
    imwrite(u(:, :, k), fullfile(outputFolder, sprintf('%s-denoised-%.3i.png', dataset, k)));
end

% Optical flow.
outputFolder = fullfile(resultfolder, removebrackets(groupname), removebrackets(dataset), 'flow');
mkdir(outputFolder);
for k=1:size(v1, 3)
    col = flowToColorV2noBoundary(cat(3, v1(:, :, k), v2(:, :, k)));
    imwrite(col, fullfile(outputFolder, sprintf('%s-flow-%.3i.png', dataset, k)));
end
end