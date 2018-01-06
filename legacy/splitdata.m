% This script takes a data folder and splits the data into smaller 
% sequences e.g. of 50 frames.
clear;
close all;
clc;

% Folder name.
datafolder = fullfile('/home/ll542/store/Dropbox (Cambridge University)/Maik and Hendrik and Carola shared/Data November 2017/');

% Output folder.
outputfolder = fullfile('/home/ll542/store/Dropbox (Cambridge University)/Maik and Hendrik and Carola shared/Data November 2017-split/');

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
        end
    end
end