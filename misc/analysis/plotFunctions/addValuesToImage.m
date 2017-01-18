function addValuesToImage(image)



for i = 1:size(image,1)
    for j = 1:size(image,2)
%     j=i;
        if image(i,j) < 0
            myCol = [.9, .4, .1];
        else
            myCol = [.1, .9, .4];
        end
        
        text(0.7+(j-1), i, sprintf('%.2f',image(i,j)), 'Color', myCol,...
            'FontSize',8);
    end
end


end
