function matrixCell = d3mat2cell(d3mat)
    % Converts a 3D matrix to a cell of 2D matrices over dimension 3
    matrixCell = squeeze(num2cell(d3mat,[1 2]))';
end