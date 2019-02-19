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

% This script creates plots for comparison of different parameters.
clear;
close all;
clc;

% Set colour-coding function.
% colourcode = @(v1, v2, len) computeColour(v1 / len, v2 / len);
% colourcode = @(v1, v2, len) flowToColorV2(cat(3, v1, v2) / len, 0);
colourcode = @(v1, v2, len) flowToColorV2(cat(3, v1, v2) / len, 10);

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

% Define folder with results.
%method = 'joint-approach';
method = 'standard-optical-flow';
resultfolder = fullfile('results', 'parameterstudies', method);

% Define figure folder.
figurefolder = fullfile('results', 'figures', 'parameterstudies', method);
mkdir(figurefolder);

fprintf('Starting analysis of folder: %s\n', datapath);
fprintf('Found %i groups.\n', length(groups));

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
        
        % Scan directory.
        folderContent = dir(fullfile(outputfolder, '*.mat'));
        fprintf('Found %i results.\n', numel(folderContent));
        if(numel(folderContent) == 0)
            continue;
        end
        
        % Load results.
        v1 = cell(numel(folderContent), 1);
        v2 = cell(numel(folderContent), 1);
        uinit = cell(numel(folderContent), 1);
        u = cell(numel(folderContent), 1);
        f = cell(numel(folderContent), 1);
        for p=1:numel(folderContent)
            R = load(fullfile(outputfolder, folderContent(p).name));
            v1{p} = R.v(:, :, :, 1);
            v2{p} = R.v(:, :, :, 2);
            uinit{p} = R.uinit;
            u{p} = R.u;
            f{p} = R.f;
            fprintf('Result %.2i', p);
            for q =1:length(R.params)
                fprintf(', param%i=%g', q, R.params{q}(p));
            end
            fprintf('\n');
        end
        
        % Create folders.
        mkdir(fullfile(figurefolder, removebrackets(groupname)));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'flow'));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'flow-first-frame'));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'flowmean'));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'denoised'));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'denoised-first-frame'));
        mkdir(fullfile(figurefolder, removebrackets(groupname), 'composition'));
        
        % Output sequence.
        f = f{1};
        [n, m, t] = size(f);
        fseq = reshape(uint8(255*f(:, :, 1:end)), n, m * t);
        %fseq = cat(3, fseq, fseq, fseq);
        imwrite(fseq, fullfile(figurefolder, removebrackets(groupname), sprintf('%s.png', dataset)));

        % Compute mean of flow.
        meanv1 = cellfun(@(x) mean(x, 3), v1, 'UniformOutput', false);
        meanv2 = cellfun(@(x) mean(x, 3), v2, 'UniformOutput', false);
        meanlen = max(cellfun(@(x, y) max(hypot(x(:), y(:))), meanv1, meanv2, 'UniformOutput', true));
        
        % Compute scaling for flow.
        len = 1;
        % len = max(cellfun(@(x, y) max(hypot(x(:), y(:))), v1, v2, 'UniformOutput', true));
        % len = max(len, meanlen);
        
        uinitf = cellfun(@(x) x(:, :, 1), uinit, 'UniformOutput', false);
        uf = cellfun(@(x) x(:, :, 1), u, 'UniformOutput', false);
        for p=1:numel(folderContent)
            v1r = reshape(v1{p}, n, m * (t - 1));
            v2r = reshape(v2{p}, n, m * (t - 1));
            col = colourcode(v1r, v2r, len);
            imwrite(col, fullfile(figurefolder, removebrackets(groupname), 'flow', sprintf('%s-flow-%.2i.png', dataset, p)));
            imwrite(colourcode(v1{p}(:, :, 1), v2{p}(:, :, 1), len), fullfile(figurefolder, removebrackets(groupname), 'flow-first-frame', sprintf('%s-flow-first-frame-%.2i.png', dataset, p)));
            
            imwrite(reshape(u{p}, n, m * t), fullfile(figurefolder, removebrackets(groupname), 'denoised', sprintf('%s-denoised-%.2i.png', dataset, p)));
            imwrite(uf{p}, fullfile(figurefolder, removebrackets(groupname), 'denoised-first-frame', sprintf('%s-denoised-first-frame-%.2i.png', dataset, p)));
            
            colmean = colourcode(meanv1{p}, meanv2{p}, len);
            imwrite(colmean, fullfile(figurefolder, removebrackets(groupname), 'flowmean', sprintf('%s-flowmean-%.2i.png', dataset, p)));
            
            uinitr = 255 * repmat(uinitf{p}, [1, 1, 3]);
            ur = 255 * repmat(uf{p}, [1, 1, 3]);
            res = cat(2, uinitr, ur, col, colmean);
            imwrite(res, fullfile(figurefolder, removebrackets(groupname), 'composition', sprintf('%s-composition-%.2i.png', dataset, p)));
        end

        % Compute flow (useful if results need to have same scaling).
        % x1 = concatenatecellarrays(v1, n, m, t - 1);
        % x2 = concatenatecellarrays(v2, n, m, t - 1);
        % meanv1 = concatenatecellarrays(meanv1, n, m, 1);
        % meanv2 = concatenatecellarrays(meanv2, n, m, 1);
        % col1 = cat(2, x1, meanv1);
        % col2 = cat(2, x2, meanv2);
        % flow = colourcode(col1, col2, len);
        % imwrite(flow, fullfile(resultfolder, removebrackets(groupname), sprintf('%s-results.png', dataset)));
        
        % Output first frame of initially denoised sequence.
        %uinit = concatenatecellarrays(uinit, n, m, 1);
        %imwrite(uinit, fullfile(resultfolder, removebrackets(groupname), sprintf('%s-denoised.png', dataset)));
        %uinit = 255 * cat(3, uinit, uinit, uinit);
        %uinit = cat(1, 255*ones(n, m, 3), uinit);

        % Output first frame of reconstructed sequence.
        %u = concatenatecellarrays(u, n, m, 1);
        %imwrite(u, fullfile(resultfolder, removebrackets(groupname), sprintf('%s-reconstructed.png', dataset)));
        %u = 255 * cat(3, u, u, u);
        %u = cat(1, 255*ones(n, m, 3), u);
    end
end

function y = concatenatecellarrays(x, n, m, t)
    y = cell2mat(cellfun(@(x) reshape(x, n, m * t), x, 'UniformOutput', false));
end