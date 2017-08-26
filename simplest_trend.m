function EsMF = simplest_trend(X)

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

% Establishing left and right averaged thresholds 
Rthbar = [2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366 2.7366];
Lthbar = [1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163 1.8163];

% ---------------------------------------------------------------
% OPTIONS.INTERP = oint;
% OPTIONS.MAXITERATIONS = 50;
imf = emd(X,'STOP',[0.2,0.5,0.05],'MAXITERATIONS',50);

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

for k = 1:L
energy_imf(k)=sum(imf(k,:).^2);
end


%% TREND ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EsMF = zeros(1,N);
for l = MInC :L
EsMF =  EsMF + imf(l,:);
end

% %% DETRENDING THE MIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% EsMF_detrended = zeros(1,N);
% for l = 1:MInC-1
% EsMF_detrended =  EsMF_detrended + imf(l,:);
% end
% 
% %% COMPUTING TEST STATISTIC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% R2 = var(X)/var(EsMF_detrended);
% 
% %% Residue
% 
% [row,column] = size(imf);
% 
% residue = imf(row,:);
% residue_index = row;

