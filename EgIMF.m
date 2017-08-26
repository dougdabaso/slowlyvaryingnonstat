function [Eg,VecEg] = EgIMF(imf)

% Author: Azadeh Moghtaderi (September 2009)
% ----------------------------------------------------------------
% This function takes the IMFs associated with a given set of data
% and calculates the energy of each IMF. It then labels those IMF 
% indexes whose associated energy have increased.
% -----------------------------------------------------------------
% INPUT: 
% imf : IMFs associated with some given data set
% -----------------------------------------------------------------
% OUTPUT:
% Eg : Energy in each IMF
% VecEg : This is a vector which marks the indexes of IMFs by their 
%         index number where the energy of that IMF is increased 
%         comparing with the energy of the previous IMF.
% -----------------------------------------------------------------

L = length(imf(:,1));
Eg = zeros(1,L) ;
for l = 1:L
    Eg(l) = sum(imf(l,:).^2) ; %Energy of the lth IMF
end


VecEg = zeros(L,1);
te = Eg(1);
for l = 2:L
    if (Eg(l) > te) 
     VecEg(l,1) = l; % mark the indexes where energy has increased
     te = Eg(l);
    else
     te = Eg(l);   
    end
end