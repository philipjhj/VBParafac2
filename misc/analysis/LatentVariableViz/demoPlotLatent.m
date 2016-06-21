%% Demostration: Plot (latent) spatial maps and associated temporal activation
% DESCRIPTIVE TEXT
clear; close all; clc
load example_data
%%
% W has dimension V x q, where each column is a spatial map
% Z has dimension q x T x subjects,
% mask is the binary mask as a vector
% mask_dim is the dimensionality of the measured 3D volume
plotComponent(W,Z,mask,mask_dim)

%% Plotting a subset or using variable arguments
subset = 1:5;
plotComponent(W(:,subset),Z(subset,:,:),mask,mask_dim,...,)
            'inConvention','Radiological',... %Input convention (i.e. data)
            'outConvention','Neurological',...%Output convention
            'TR',2.49,... %Time resolution in secounds
            'FontSize',18,... % Self-explanatory. Note titles has Fontsize+2 
            'LineWidth',2,... % Width of the temporal activation
            'save','sdas',... %Saves .png and .eps to the given directory, if '' figures aren't saved
            'Position',[0 0 800 425]); %Figure position lower left and upper right cornor
