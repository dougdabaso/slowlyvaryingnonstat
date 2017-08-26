%% Stationarity test based on trend analysis of the time-marginal
%% distributoon

% Input:  
%         x - signal 
%         N_btstrp - Number of bootstraps used in the analysis 
%         Nh - Size of the hermite window
%         n_tapers - number of tapers       
%
% Output:
%        
% res_test - stationarity test result (1 - nonstat, 0 - stat)
% R2_signal - trend importance estimation of the signal itself
% gev_thresh - GEV threshold 
% gev_hat - vector containing the GEV distribution parameters
% R2_bootstrap - distribution of the trend importance estimator in the stationary data
% B2 - Zempleni statistic
% TI_lb_signal - index of nonstat

% chosen_window:
%
% 1) Hermite tapers
% obs: default number of tapers: 10
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


function [res_test,R2_signal,gev_thresh,gev_hat,R2_bootstrap,B2,TI_lb_signal,Eg_vec] = stat_test_TI_analysis2(x,N_btstrp,Nh,chosen_window,fa_rate_array);


% if Nh=1 we don't need windows
if Nh == 1
chosen_window = 0;
end

% 1) Setting  variables and parameters that will be used in this function
%fa_rate = 0.05;

% 2) Separate the blocks to perform the block bootstrapping and generate
% virtual realizations of the signal. Possible correlations in the series
% are preserved within a block

% nonoverlapping blocks, setting the blocks to be scrambled
desired_length = length(x);
bl = ceil((desired_length)^(0.25));
number_blocks = floor(desired_length/bl);
[y]=segmcirc(x,bl,bl,number_blocks);
[rows,cols] = size(y); % blocks in the variable y

  
% 3) Compute the trend importance in the original signal

% 3.1) Computing the time marginal
time_marginal = mostrecent_marginal_in_time(x,Nh,chosen_window);
[EsMF,Eg,MInC] = simpler_trend(time_marginal);
R2_signal = var(time_marginal)/var(time_marginal-EsMF);

% 3.3) Computing the index of nonstationarity
TI_lb_signal = TI_index_of_nonstat(Eg,MInC);
%TI_lb_signal = 0;

Eg_vec = Eg;

clear time_marginal
clear Eg
clear MInC

% 4) Generate resamples and compute the test statistic  
for i = (1:N_btstrp)   
r = ceil(cols.*rand(cols,1));
resampled_x = y(:,r);
resampled_x = resampled_x(:)'; 
vector_resampled_x(i,:) = resampled_x;
end

% 4.1) Computing the time marginal from the resamples 
%parfor i_in = (1:N_btstrp) 
for i_in = (1:N_btstrp) 
time_marginal = mostrecent_marginal_in_time(vector_resampled_x(i_in,:),Nh,chosen_window);
R2_bootstrap(i_in) = var(time_marginal)/var(time_marginal-simplest_trend(time_marginal));
end


% 6) Performing the hypothesis test
for i_alpha = 1:length(fa_rate_array)

fa_rate = fa_rate_array(i_alpha);
[gev_hat,parmci] =  gevfit(R2_bootstrap);
gev_thresh=gevinv(1-fa_rate,gev_hat(1),gev_hat(2),gev_hat(3));
res_test.(sprintf('i_alpha_%d',i_alpha)) = (R2_signal>gev_thresh);

end

% 7) Computing the Zempleni statistic
dist = R2_bootstrap(:);
dist = sort(dist);
z = gevcdf(dist,gev_hat(1),gev_hat(2),gev_hat(3));    
t = 1:N_btstrp;
S1 = sum(((2*t)-1)*log(1-z(N_btstrp+1-t))/N_btstrp);
S2 = sum(z/N_btstrp);
B2 = (N_btstrp/2)-S1-S2;


end


