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
%addpath(genpath(cd))

%folder name
mainFolder = fullfile('/home/ll542/store/Dropbox (Cambridge University)/Maik and Hendrik and Carola shared/Data November 2017 Reproduce/EB1 Data', 'control');
%mainFolder = [datapath, '01_control'];
%mainFolder = 'data';

% Add all subfolders.
y = dir(mainFolder);
y = y(~cellfun(@(x) strcmp(x, '.git') || strcmp(x, '.') || strcmp(x, '..'), {y.name}));
listFolders = y([y.isdir]);

for i=1:numel(listFolders)
    folder = [mainFolder,filesep,listFolders(i).name,filesep];
    result = dataPipelineFunction_flow(folder);
end