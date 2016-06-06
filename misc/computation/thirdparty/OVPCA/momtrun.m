function [mv,var,nk]=tN(m,s,a,b);
% function evaluates moments of truncated normal distribution

alf	= (a-m)./(sqrt(2)*s);
bet	= (b-m)./(sqrt(2)*s);

pom	= (erf(bet) - erf(alf))* sqrt(pi/2);
gam	= (exp(-bet.^2) - exp(-alf.^2)) ./ pom;
del	= (b.*exp(-bet.^2) - a.*exp(-alf.^2)) ./ pom;

mv	= m - s.*gam;
var	= s.^2 + m.*mv - s.*del;
nk	= s .* pom;

if any(mv<a)
	disp('momtrun:');keyboard
end
