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

% This script takes a data folder and splits the data into smaller 
% sequences e.g. of 50 frames.
% This is required for the joint approach due to memory issues.
clear;
close all;
clc;

% Folder name.
datafolder = fullfile('/home/ll542/Dropbox (Cambridge University)/Maik and Hendrik and Carola shared/Data November 2017/');

% Output folder.
outputfolder = fullfile('/home/ll542/Dropbox (Cambridge University)/Maik and Hendrik and Carola shared/Data November 2017-split/');

% Sequence length.
maxlength = 50;

% Scan subfolders.
y = dir(datafolder);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groupfolders = y([y.isdir]);

fprintf('Scanning data folder: %s\n', datafolder);

% Run through groups.
for k=1:length(groupfolders)
    fprintf('Scanning group: %s\n', groupfolders(k).name);
    groupfolder = fullfile(datafolder, groupfolders(k).name);
    
    % Scan datasets.
    y = dir(groupfolder);
    y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
    datasets = y([y.isdir]);
    
    % Run through datasets.
    for l=1:length(datasets)
        fprintf('Scanning dataset: %s\n', datasets(l).name);
        datasetfolder = fullfile(datafolder, groupfolders(k).name, datasets(l).name);
        
        % Scan and sort images.
        files = dir(fullfile(datasetfolder, '*.tif'));
        files = sort({files.name});
        
        % Determine number of sequences.
        n = ceil(length(files) / maxlength);
       
        % Create folder and copy files.
        for s=1:n
            sequencename = sprintf('%s_%i', datasets(l).name, s);
            targetfolder = fullfile(outputfolder, groupfolders(k).name, sequencename);
            mkdir(targetfolder);
            
            % Copy files.
            offset = maxlength*(s-1);
            for f=1:maxlength
                copyfile(fullfile(datasetfolder, files{f+offset}), fullfile(targetfolder, files{f+offset}));
            end
            
            % Copy segmentation.
            mkdir(fullfile(targetfolder, 'images'));
            copyfile(fullfile(datasetfolder, 'images', 'segmentationMap.png'), fullfile(targetfolder, 'images', 'segmentationMap.png'));
            
            % Copy region maps if present.
            if(isfile(fullfile(datasetfolder, 'images', 'segmentationMap1.png')) && isfile(fullfile(datasetfolder, 'images', 'segmentationMap2.png')))
                copyfile(fullfile(datasetfolder, 'images', 'segmentationMap1.png'), fullfile(targetfolder, 'images', 'segmentationMap1.png'));
                copyfile(fullfile(datasetfolder, 'images', 'segmentationMap2.png'), fullfile(targetfolder, 'images', 'segmentationMap2.png'));
            end
        end
    end
end