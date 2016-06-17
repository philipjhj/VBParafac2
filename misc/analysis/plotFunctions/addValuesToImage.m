


function addValuesToImage(image)



for i = 1:size(image,1)
    for j = 1:size(image,2)
        text(0.8+(j-1), i, sprintf('%.2f',image(i,j)), 'Color', [.9, .4, .1],...
            'FontSize',8);
    end
end


end
