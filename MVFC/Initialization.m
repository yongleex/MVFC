%% Initialization 
% clear the work space and close all figures 
clear; close all;clc;

% Add all the folders into Matlab path
currentFolder = pwd;
addpath(genpath(currentFolder));

% Change the plot model in case of drawing crash
opengl('save', 'software');

% Finished
fprintf('Initialized!\n')
