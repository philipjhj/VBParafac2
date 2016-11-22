function score = congruenceScore(x,y)
    xy = x'*y;
    xxsq = sqrt(x'*x);
    yysq = sqrt(y'*y);
    
    score = xy/(xxsq*yysq);
    
%     xxsq = sqrt(diag(xy));
%     xxyy = xxsq*xxsq';
%     
%     if numel(xy)>1
%         score = triu(xy./(xxyy),1);
%     else
%         score = xy./(xxyy);
%     end
end