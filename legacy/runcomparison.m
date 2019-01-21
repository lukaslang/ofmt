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

% This script runs a comparison between optical flow implementations.
clear;
close all;
clc;

% Load input sequence.
dataset = 'RubberWhale';
forig(:, :, 1) = im2double(imread(fullfile('benchmark', 'other-data-gray', dataset, 'frame10.png')));
forig(:, :, 2) = im2double(imread(fullfile('benchmark', 'other-data-gray', dataset, 'frame11.png')));

% Scale images.
scale = 0.5;
for k=1:size(forig, 3)
    f(:, :, k) = imresize(forig(:, :, k), scale, 'bilinear', 'Antialiasing', false);
end
[n, m, t] = size(f);

% Set output folder.
outputfolder = fullfile('results', 'figures', 'comparison');
mkdir(outputfolder);

%% Joint model.

% Set regularisation parameter.
alpha = 1e-12;
beta = 0.00001;
gamma = 0.001;

tic;
[uinit, u, v] = runjointmodel(f, alpha, beta, gamma);
toc;
v = v(:, :, 1:end-1, :);

figure;
img = flowToColorV2(squeeze(v));
imagesc(img);
axis image;
drawnow();

figure;
imagesc(abs(uinit(:, :, 1) - f(:, :, 1)));
axis image;
colorbar;
title('uinit - f');
drawnow();

figure;
imagesc(abs(u(:, :, 1) - f(:, :, 1)));
axis image;
colorbar;
title('u - f');
drawnow();

%% Standard optical flow.

% Set regularisation parameter.
alpha = 0.0005;
beta = 0.005;

tic;
v = runstandardof(f, alpha, beta);
toc;

figure;
img = flowToColorV2(squeeze(v));
imagesc(img);
axis image;
drawnow();

%% Functions.

function [uinit, u, v] = runjointmodel(f, alpha, beta, gamma)
    temporalSmoothness = 0;
    mainJoint = jointModelLargeScale(f, alpha, beta, gamma, 'temporalSmoothness', temporalSmoothness, 'verbose', 1);
    mainJoint.doWarping = false;
    mainJoint.opticalFlowTerm = 'classic';
    mainJoint.numMinMainIt = 10;
    mainJoint.numMainIt = 10;
    
    % Init and get denoised sequence.
    mainJoint.init; 
    uinit = mainJoint.u;

    % Run algorithm.
    mainJoint.run;

    % Recover solution.
    u = mainJoint.u;
    v = mainJoint.v;
end

function v = runstandardof(f, alpha, beta)
    [n, m, t] = size(f);
    
    % Define optical flow problem.
    p = @vecof3dl2tv;

    % Set algorithm parameters.
    tau = 1/sqrt(8);
    sigma = 1/sqrt(8);

    % Define termination criterion.
    term = @(iter, p, pprev, tau, sigma) pdresidual(p, pprev, tau, sigma) < 1e-6;

    % Define verbosity, logging, set plotting callback.
    alg = pdhg(tau, sigma, 1000, term, 1000, @logenergy);

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
    
    % Set up problem.
    p = p(fx, fy, ft, alpha, beta, Dx, Dy, Dt, v, y);

    % Run algorithm.
    stats = alg.run(p);

    % Recover solution.
    [u, v] = p.solution;
    v = cat(4, u, v);
end

function logfun(iter, p, pprev, tau, sigma, plot)
    if(plot)
        logplotenergy(iter, p, pprev, tau, sigma);
    else
        logenergy(iter, p, pprev, tau, sigma);
    end
end