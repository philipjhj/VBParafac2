function [avg_measure, idx, matched_measure,sorted_measures,idx_rows_final] = greedy_component_match(measures)
% Match
% measures: is a D x D matrix, where the measure is between the estimated
% components (rows) and true components (columns)
measures_temp=measures;
[D1,D2] = size(measures);

% assert(all(size(measures)==D),'The input was not a square matrix!')
assert(all(measures(:)>0),'A measure was negative, this is not supported.')

idx = nan(D2,1);
matched_measure = nan(D2,1);
for d2 = 1:D2
    % Find the value and column index of the maximum value
    [max_meas,tmatch] = max(max(measures,[],1));
    % Find the estimated components which produced this match
    [~,test] = max(measures(:,tmatch)); 
    
    % save it
    idx(d2) = tmatch;
    matched_measure(d2) = max_meas;
    
    % Set this component as matched
    measures(:,tmatch) = 0;
%     measures(test,:) = 0;
end

%Sort rows
measures_temp=measures_temp(:,idx);
idx_rows_final=1:D1;
for d1 = 1:min(D1-1,D2)
    [~,idx_rows]=sort(measures_temp(d1:end,d1),'descend');
    measures_temp(d1:end,:)=measures_temp(idx_rows+d1-1,:);
    idx_rows_final(d1:end)=idx_rows_final(idx_rows+d1-1);
end
sorted_measures=measures_temp;

avg_measure = mean(matched_measure);

end