function rgb = HInt2RGB (Hue, Int)

% This function determines the rgb output corresponding to input values Hue
% and Int. Hue has to be an integer 1 - 12 and Int has to be an integer
% comprised between 1 and 100 .
% Hue is the Hue according to the HSV convention
% (https://en.wikipedia.org/wiki/HSL_and_HSV) and indicates the "type of
% color" in 30 deg increment following the HSV map (i.e. 1 corresponds to
% 0 deg and is "red", 2 corresponds to 30 deg and is "yellow-red", 3
% corresponds to 60 deg and is "yellow".  More "different" colors are
% located at highest angular distances, i.e. the highest contrast to the
% red color (0 deg) is the cyan (180 deg).
% The value Int (integer between 1 and 100) is a measure of darkness 
% of the color when exported to grayscale. The values have been 
% calculated trough minimization for all the different swatches (hue 
% colors), to match the intensity of pure red (HSV=[0 1 1]). 
% Note that pure red is not as black as black (i.e. rgb2gray([1 0 0])=0.2989),
% so values with Int=100 levels will look gray. Perfect for contrast with 
% actual black lines 

load ('HSVBWColormap.mat');

rgb=hsv2rgb([(Hue-1)/12,Ss(Int,Hue),Vs(Int,Hue)]);
end