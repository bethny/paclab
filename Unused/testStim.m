%% ADD CODE BASES TO PATH
addpath(genpath('~/code/pac/paclab'));

%%
clear all

lambda = 20; % density of the gabor stripes
thetaStep = [-45,-30,-20,-10,0,10,20,30,45];
imSize = 100;

makeGabor(lambda,imSize,thetaStep)