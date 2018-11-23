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
resultfolder = fullfile('results', 'parameterstudies');
mkdir(resultfolder);

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
            R = load(fullfile(outputfolder, folderContent(p).name), 'f', 'uinit', 'u', 'v');
            v1{p} = R.v(:, :, :, 1);
            v2{p} = R.v(:, :, :, 2);
            uinit{p} = R.uinit;
            u{p} = R.u;
            f{p} = R.f;
        end

        % Plot results.
        f = f{1};
        [n, m, t] = size(f);
        x1 = concatenatecellarrays(v1, n, m, t);
        x2 = concatenatecellarrays(v2, n, m, t);
        fseq = reshape(uint8(255*f(:, :, 1:end)), n, m * t);
        fseq = cat(3, fseq, fseq, fseq);
        
        uinit = concatenatecellarrays(uinit, n, m, t);
        uinit = 255 * cat(3, uinit, uinit, uinit);
        uinit = cat(1, fseq, uinit);
        
        u = concatenatecellarrays(u, n, m, t);
        u = 255 * cat(3, u, u, u);
        u = cat(1, fseq, u);
        
        flow = cat(1, fseq, flowToColorV2noBoundary(cat(3, x1, x2)));
        img = cat(2, uinit, u, flow);
        imwrite(img, fullfile(resultfolder, removebrackets(groupname), sprintf('%s-results.png', dataset)));
        figure;
        imagesc(img);
        axis image;
        title('Results.');
    end
end

function y = concatenatecellarrays(x, n, m, t)
    y = cell2mat(cellfun(@(x) reshape(x, n, m * t), x, 'UniformOutput', false));
end