% Note:
% (1) adjust the folder variable
% (2) create a subfolder in folder named as 'images'
% (3) take one input image and cut out the region of interest, et the outer
% region to 0 and save the image as 'segmentationMap.png' in the subfolder
% 'images'
% (4) run this script
%
clear;
close all;
clc;

% Set output folder.
resultfolder = 'results';

% Add all subfolders.
y = dir(datapath);
y = y(~cellfun(@(x) strcmp(x, '.') || strcmp(x, '..'), {y.name}));
groups = y([y.isdir]);

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
        fprintf('Dataset: %s\n', fullfile(groupname, dataset));
        
        % Create output folder.
        outputfolder = fullfile(resultfolder, groupname, dataset);
        mkdir(outputfolder);
        
        % Create plots.
        tic;
        dataPipelineFunction_plots(datafolder, outputfolder);
        toc;
    end
end