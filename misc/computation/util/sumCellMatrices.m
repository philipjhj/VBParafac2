function sumMat = sumCellMatrices(myCell)
    % Sums matrices stored in a cell array
    sumMat = sum(cat(3,myCell{:}),3);
end