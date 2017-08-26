function [Rt,VecRt] = RzcnIMF(imf,Lthbar,Rthbar)

% Author: Azadeh Moghtaderi (September 2009)
% -----------------------------------------------------------------
% This function uses the zero crossing numbers (ZCN) of IMFs 
% associated with a given data set and calculates the ratio between 
% every two consecutive ZCN. It then labels by one those IMF indexes 
% whose associated ratio is outside of the given threshhold.
% -----------------------------------------------------------------
% INPUT: 
% imf : IMFs associated with some given data set
% Lthbar : Averaged left threshold.
% Rthbar : Averaged right threshold.
% -----------------------------------------------------------------
% OUTPUT:
% Rt : The vector of ratios.
% VecRt : This is a vector which marks the indexes of IMFs by one 
%         where the ratio between two consecutive ZCN is outside of 
%         the given threshold.
% -----------------------------------------------------------------

L = length(imf(:,1));
N = length(imf(1,:));
indz = zeros(1,L);
t = 1:N;
for l = 1:L

   % I'm computing the maxima, minima and zero crossing indices for the
   % first imf given in imf(1,:)
   [indmin, indmax, indzer] = extr(imf(l,:),t);
   indz(l) = length(indzer);
end


% If the last IMF has ZCN =0, we add a little bit to it not to have problem
% when dividing.
%indzt = zeros(1,L);
if (indz(L) == 0)
   indzt(1:L-1) =  indz(1:L-1); 
   indzt(L) = indz(L)+eps;
else 
   indzt(1:L) =  indz(1:L);  
end


% Calculating the ratios between ZCN of every two consecutive IMF
Rt = zeros(1,L-1);
for j = 1:L-1
    Rt(j) = indzt(j) ./ indzt(j+1);
end

% ----------------------------------------------------------
VecRt= zeros(L,1);
if L-1 < 14
    for l = 1:L-1
    if (Rt(l) > Rthbar(l))  || (Rt(l) < Lthbar(l))
        VecRt(l+1,1) = l+1;
    end
    end
else 
    for m = 1:14
        if (Rt(m) > Rthbar(m))  || (Rt(m) < Lthbar(m))
            VecRt(m+1,1) = m+1;
        end
    end
    for m = 15:L-1
        if (Rt(m) > Rthbar(14))  || (Rt(m) < Lthbar(14))
            VecRt(m+1,1) = m+1;
        end
    end
end



