% Fast estimation of the time marginal

% input:
%
% x: original signal
% Nh: number of points of the chosen window
%
% chosen_window:
%
% 0) No window
% obs: time marginal = abs(x)^2
%
% 1) Hermite tapers
% obs: default numbrer of tapers: 10
%
% 2) Hanning window
% obs: symmetric Hanning window
%
% 3) Kaiser window
% obs: free parameter (beta) set to 0.5
%
% 4) Hamming window
% obs: symmetric Hamming window
%
% 5) Bartlett window
%
% 6) triangular window
%
% 7) rectangular window
%
% 8) Gaussian window
% obs: alpha = 2.5


function time_marginal_fast = mostrecent_marginal_in_time(x,Nh,chosen_window);


if chosen_window == 0
time_marginal_fast=x.*x;
    
elseif chosen_window == 1
    
% Hermite window
[h,Dh,tt] = hermf(Nh,10,6);
for i = 1:size(h,1)
conv_with_taper(i,:)=(conv(x',h(i,:))).^2;
end
time_marginal_fast = sum(conv_with_taper);
    

elseif chosen_window == 2

% Hanning window
time_marginal_fast=(conv(x',hanning(Nh))).^2;
    
    
elseif chosen_window == 3
    
% Kaiser window    
time_marginal_fast=(conv(x',kaiser(Nh))).^2;
    
elseif chosen_window == 4

 % Hamming window  
 time_marginal_fast=(conv(x',hamming(Nh))).^2;
    
elseif chosen_window == 5
    
 % Bartlett window  
 time_marginal_fast=(conv(x',bartlett(Nh))).^2;    
    
elseif chosen_window == 6
    
 % triangular window  
 time_marginal_fast=(conv(x',triang(Nh))).^2;
 
elseif chosen_window == 7
    
% rectangular window  
 time_marginal_fast=(conv(x',rectwin(Nh))).^2;    
 
elseif chosen_window == 8
    
 % Gaussian window  
 time_marginal_fast=(conv(x',gausswin(Nh))).^2;    
    
end

if size(time_marginal_fast,1) > 1
    time_marginal_fast = time_marginal_fast';
end

