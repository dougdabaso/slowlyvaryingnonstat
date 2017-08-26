function [CompMat,InC,MInC,EsMF,Rt,Eg,L,imf,Lthbar,Rthbar] = EgRtComp(X,oint)

% Author: Azadeh Moghtaderi (September 2009)
% -----------------------------------------------------------------
% This function compares the result from EgIMf.m and RzcnIMF.m to 
% make a conclution about how many modes to select for trend estimate.
% -----------------------------------------------------------------
% INPUT: 
% X : data set of interest
% oint : Interpolation option in EMD chosen to be either 'linear' or 'spline'
% -----------------------------------------------------------------
% OUTPUT:
% CompMat: A matrix which shows the outcomes of EgIMf.m and RZCNIMF.m
%          together.      
% MInC :  Minimum index of nonzero entry common in both EgIMf.m and 
%         RzcnIMF.m
% EsLMF : Estimated modulating (trend) function.
% ------------------------------------------------------------------------

% loading left and right averaged thresholds 
load Lthbar-spline-21ex-a0.09-flat.mat
load Rthbar-spline-21ex-a0.09-flat.mat

% ---------------------------------------------------------------
OPTIONS.INTERP = oint;
OPTIONS.MAXITERATIONS = 50;
imf = emd(X,OPTIONS);

% ----------------------------------------
[Eg,VecEg] = EgIMF(imf);
[Rt,VecRt] = RzcnIMF(imf,Lthbar,Rthbar);
L = length(imf(:,1));
N = length(imf(1,:));

% Create one matrix which compares the energy and ratio result
CompMat = zeros(L,2);
CompMat(:,1) = VecEg;
CompMat(:,2) = VecRt;

% ----------------------------------------
k = 1;
InC=[];
for l = 1:L
   if (CompMat(l,1) ~= 0) && (CompMat(l,2) ~= 0)
     InC(k) = l; % indexes with nonzero entries in both energy and ratio
     k = k+1;
   end
end

if ~isempty(InC)
    MInC  = min(InC);
else
    MInC = length(CompMat);
end

% ----------------------------------------------------------

% estimation the trend 
EsMF = zeros(1,N);

for l = MInC :L
EsMF =  EsMF + imf(l,:);
end

