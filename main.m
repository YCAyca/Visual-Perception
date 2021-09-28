clear all;
close all;
clc;

format short

%% DLT

% Loading 3D World Points and 2D projected points (Ground Truth Points for projection)coming from previous lab
World_PTS = load('pts3D.txt'); 
GT_PTS = load('PTS_.txt');
total_point_count = 6; % We will use 6 points to do the DLT Calibrations

DLT(World_PTS, GT_PTS,total_point_count);

