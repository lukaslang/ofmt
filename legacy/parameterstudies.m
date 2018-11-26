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

% This script runs parameter studies for the joint approach.
clear;
close all;
clc;

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Define and create folder with results.
resultfolder = fullfile('results', 'parameterstudies', 'joint-approach');
mkdir(resultfolder);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Set parameter range (alpha, beta, gamma) and generate combinations.
params = {[5e-3, 1e-2, 5e-2], [1e-2, 5e-2, 1e-1], [1e0, 1e1, 1e2]};
[param1, param2, param3] = ndgrid(params{:});
ncombs = numel(param1);

% Run through all groups.
for k=1:length(groups)
    groupname = groups(k).name;
    % Run through all datasets.
    y = dir(fullfile(datapath, groupname));
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    for l=1:length(datasets)
        dataset = datasets(l).name;
        datafolder = fullfile(datapath, groupname, dataset, filesep);
        outputfolder = fullfile(resultfolder, removebrackets(groupname), dataset);
        
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        mkdir(outputfolder);
        
        % Load dataset.
        folderContent = dir(fullfile(datafolder, '*.tif'));
        for q=1:numel(folderContent)
            f(:, :, q) = im2double(imread(fullfile(datafolder, folderContent(q).name)));
        end
        
        % Compute and save results.
        for p=1:ncombs
            alpha = param1(p);
            beta = param2(p);
            gamma = param3(p);
            fprintf('Parameter combination %.2i/%.2i.\n', p, ncombs);
            tic;
            [uinit, u, v] = runjointmodel(f, alpha, beta, gamma);
            toc;
            u = u(:, :, 1:end-1);
            v = v(:, :, 1:end-1);
            save(fullfile(outputfolder, sprintf('results-%.2i.mat', p)), 'f', 'uinit', 'u', 'v', 'alpha', 'beta', 'gamma');
        end
    end
end

% Define and create folder with results.
resultfolder = fullfile('results', 'parameterstudies', 'standard-optical-flow');
mkdir(resultfolder);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

% Set parameter range (alpha, beta, gamma) and generate combinations.
params = {[5e-3, 1e-2], [1e-2, 1e-1], [1e-4, 1e-3], [1e-3, 1e-2]};
[param1, param2, param3, param4] = ndgrid(params{:});
ncombs = numel(param1);
    
% Run through all groups.
for k=1:length(groups)
    groupname = groups(k).name;
    % Run through all datasets.
    y = dir(fullfile(datapath, groupname));
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    for l=1:length(datasets)
        dataset = datasets(l).name;
        datafolder = fullfile(datapath, groupname, dataset, filesep);
        outputfolder = fullfile(resultfolder, removebrackets(groupname), dataset);
        
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        mkdir(outputfolder);
        
        % Load dataset.
        folderContent = dir(fullfile(datafolder, '*.tif'));
        for q=1:numel(folderContent)
            f(:, :, q) = im2double(imread(fullfile(datafolder, folderContent(q).name)));
        end
        
        % Compute and save results.
        for p=1:ncombs
            alpha1 = param1(p);
            beta1 = param2(p);
            alpha2 = param3(p);
            beta2 = param4(p);
            fprintf('Parameter combination %.2i/%.2i.\n', p, ncombs);
            tic;
            [uinit, u, v] = runstandardof(f, alpha1, beta1, alpha2, beta2);
            toc;
            save(fullfile(outputfolder, sprintf('results-%.2i.mat', p)), 'f', 'uinit', 'u', 'v', 'alpha1', 'beta1', 'alpha2', 'beta2');
        end
    end
end

function [uinit, u, v] = runjointmodel(f, alpha, beta, gamma)
    temporalSmoothness = 0.5;
    mainJoint = jointModelLargeScale(f, alpha, beta, gamma, 'temporalSmoothness', temporalSmoothness, 'verbose', 1);
    mainJoint.doWarping = false;
    mainJoint.opticalFlowTerm = 'classic';
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


function [f, u, v] = runstandardof(fdelta, alpha1, beta1, alpha2, beta2)
    [n, m, t] = size(fdelta);

    % Define denoising problem.
    p = @denoise3dl2tv;

    % Set algorithm parameters.
    tau = 1/sqrt(8);
    sigma = 1/sqrt(8);

    % Define termination criterion.
    term = @(iter, p, pprev, tau, sigma) pdresidual(p, pprev, tau, sigma) < 1e-6;

    % Define verbosity, logging, set plotting callback.
    alg = pdhg(tau, sigma, 50, term, 100, @logenergy);

    % Set up problem.
    p = p(fdelta, alpha1, beta1);

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
    
    % Define optical flow problem.
    p = @vecof3dl2tv;

    % Set algorithm parameters.
    tau = 1/sqrt(8);
    sigma = 1/sqrt(8);

    % Define termination criterion.
    term = @(iter, p, pprev, tau, sigma) pdresidual(p, pprev, tau, sigma) < 1e-6;

    % Define verbosity, logging, set plotting callback.
    alg = pdhg(tau, sigma, 50, term, 100, @logenergy);

    % Compute partial derivatives.
    fvec = f(:);
    [Dx, Dy, ~] = vecderiv3dc(m, n, t-1, 1, 1, 1);
    fx = reshape(Dx * fvec(1:n*m*(t-1)), n, m, t-1);
    fy = reshape(Dy * fvec(1:n*m*(t-1)), n, m, t-1);
    [~, ~, Dt] = vecderiv3dfw(m, n, t, 1, 1, 1);
    ft = reshape(Dt * fvec, n, m, t);
    ft = ft(:, :, 1:end-1);

    % Set up problem.
    p = p(fx, fy, ft, alpha2, beta2);

    % Run algorithm.
    stats = alg.run(p);

    % Recover solution.
    [u, v] = p.solution;
end