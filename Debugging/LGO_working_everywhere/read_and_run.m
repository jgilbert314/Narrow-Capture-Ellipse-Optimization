
clear all;
close all;
clc;

format long;

%%
% ==============================

% Specify N traps:

NN = 8;

% Set generic angles:

theta = pi/2*ones(4,1);
phi = pi/8*[1:8]';

Ntraps=size(theta,1);

proc_MFPT = 'eval_H_MFPT';
new_angles_MFPT = optimize_locations([theta;phi],2,proc_MFPT)
h1 = plot_traps_sph(new_angles_MFPT,0, 0, 'blue', 20);
dist_M = sorted_pairwise_dist(new_angles_MFPT);

% %%
% 
% proc_Coulomb='eval_H_Coulomb';
% new_angles_COU = optimize_locations([theta;phi],2,proc_Coulomb)
% h2= plot_traps_sph(new_angles_COU,0, 0, 'green', 20);
% dist_C=sorted_pairwise_dist(new_angles_COU);
%  

%dist_M-dist_C

