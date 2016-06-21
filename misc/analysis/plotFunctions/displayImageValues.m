function displayImageValues(image,titleText,colorInterval)
% displayImageValues: displays image with its values annotated
% Input:
%       Required:
%           image: image to be displayed
%           
%       Optional:
%           titelText: title string of image
%           colorInterval: [Min Max] values to use on colorbar
%
% If the image values are to be compared with other images, colorInterval
% can be used to fix the colorscale. If not specified, the default
% colorscale of colorbar are used.
%
% Philip J. H. Jørgensen, June 2016

imagesc(image)

if ~isempty(colorInterval)
caxis(colorInterval); 
end

% colorbar
title(titleText)
addValuesToImage(image)
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
end