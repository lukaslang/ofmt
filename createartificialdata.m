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
% This script creates noisy test data and writes it to 'data/artificial/01_000'.
clear;
close all;
clc;

% Set number of dots.
ndots = 20;

% Set number of time steps.
t = 10;

% Set sigma.
sigma = 3;

% Set subsampling.
h = 0.1;

% Set image size.
n = 50;
m = 50;

% Create initial positions of points.
X = m * rand(ndots, 1) / h;
Y = n * rand(ndots, 1) / h;

% Create positions of randomly moving points.
%X = cumsum([X, sigma * zeros(ndots, t-1)], 2);
Y = cumsum([Y, sigma * ones(ndots, t-1)], 2);

X = cumsum([X, sigma * randn(ndots, t-1)], 2);
%Y = cumsum([Y, sigma * randn(ndots, t-1)], 2);

% Round positions.
X = round(X);
Y = round(Y);

% Sample brightness uniformly.
b = rand(ndots, 1);

% Create image.
f = zeros(n / h, m / h, t);
for frame=1:t
    for k=1:ndots
        if(1 <= Y(k, frame) && Y(k, frame) <= n / h && 1 <= X(k, frame) && X(k, frame) <= m / h)
            f(Y(k, frame), X(k, frame), frame) = b(k);
        end
    end
    f(:, :, frame) = imgaussfilt(f(:, :, frame), 12, 'Padding', 'symmetric');
end

% Resize and rescale image.
fdelta = imresize3(f, [n, m, t]);
fdelta = fdelta / max(fdelta(:));

% Add noise and rescale.
fdelta = imnoise(fdelta, 'gaussian', 0, 0.05);
fdelta = fdelta / max(fdelta(:));

% Write dataset.
outputfolder = fullfile('data', 'artificial', '01_000');
mkdir(outputfolder);
for frame=1:t
	imwrite(uint8(255*fdelta(:, :, frame)), fullfile(outputfolder, sprintf('image-%.3i.tif', frame)), 'tif');
end

% Create segmentation map.
seg = ones(n, m, 'uint8');
outputfolder = fullfile(outputfolder, 'images');
mkdir(outputfolder);
imwrite(seg, fullfile(outputfolder, 'segmentationMap.png'), 'png');