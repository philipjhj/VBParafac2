function [Ka,ml,Kx,om,hist]=opca_fast(lD,p,N,approx,vr,ir);
% function [mA,mL,mX,ms,hist]=opca_fast(D,);
%
% function returns variational values of decomposed components
% using approximations
%
% approx == 1       Mardia 1977 (default)
%           2       Bessel/Bessel
%           3       Woods 2002

% start with svd
if nargin<4
	approx	= 5;
end
if size(approx,2)<2, approx=approx*[1 1]; end
% if no vector of ranks is specified, then consider it to be full
if nargin<6
	ir	= 0;
end
if nargin<5
	vr	= 1:p-1;
end

if nargout>4,	HIST	= 1; else HIST	= 0; end;

% scaling
l0	= lD/sqrt(lD'*lD);
lDorig	= lD;

for r = vr
	% Step 1, ML estimate
	lD	= l0(1:r);
	Ka	= ones(r,1);
	Kx	= ones(r,1);
	ml	= lD;

	if ir>0
		l0i	= l0(ir:end);
	else
		l0i	= l0(r+1:end);
	end
	om0	= N*(p-r)/(l0i'*l0i);
	lu	= sqrt(1./(1:r)');
	om	= om0;

	if HIST
		hist.ml	= l0(1:r);
	end
	hist.om	= om0;

	% Step 2, iterate
	i	= 2;
	diffom	= 1;
	while (diffom~=0 & i<500)	% converge really fast, ~= is affordable

		Ka	= hyperg(p/2, om*Kx.*lD.*ml,approx(1));
		Kx	= hyperg(N/2, om*Ka.*lD.*ml,approx(2));
		m	= Ka.*Kx.*lD;
		[ml,mvar]	= momtrun(m,sqrt(1/om),  0,lu);
		om	= (p*N)/(l0'*l0 - 2*sum(Kx.*ml.*Ka.*lD) + sum(mvar));

		if HIST
			hist.ml(:,i)	= ml;
		end
		hist.om(i)	= om;
		diffom	= om-hist.om(i-1);
		i	= i+1;
	end
	r

	% save values for further use
	B.Ka{r}	= Ka;
	B.Kx{r}	= Kx;
	B.ml{r} = ml;
	B.mvar{r}	= mvar;
	B.om{r}	= om;
	B.it{r}	= i;
	B.lu{r}	= lu;
	% Bayes rank sel
	if length(vr)>1
		% additional info
		% values of: Ka, Kx, ml, om, are assumed to be fixes
		[poma,VKa,lFa]	= hyperg(p/2, om*Kx.*lD.*ml,3);
		[pomx,VKx,lFx]	= hyperg(N/2, om*Ka.*lD.*ml,3);
		[poml,pomvar,nk]	= momtrun(m,sqrt(1/om),  0,lu);
		b	= N*p/2/om;

		lnB(r)	= (-r/2*log(pi)+gammaln(r)+gammaln(r/2+1)+r*log(2)...
				  +lFa +lFx -2*om*sum(lD.*Ka.*Kx.*ml)...
				  +sum(log(nk)) + 1/2*(m'*m - 2*ml'*m + sum(mvar))...
				  -(p*N/2)*log(b));
		if HIST
			B.VKa{r}	= VKa;
			B.VKx{r}	= VKx;
		end
	end

end % of ranks

if length(vr)>1
	% rank selection
	exB	= exp(lnB - max(lnB));
	fr	= exB/sum(exB);

	hist.fr	= fr;
	% results
	[pr,r]	= max(fr);
end
