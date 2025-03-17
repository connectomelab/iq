%% Kuramoto test------after obtained the optimal pair of K and V, then use this pair of K and V values run the kuramoto model to get the metastability of 280 subjects
clear;
%rng(1);

% compile the mex code:
%mex -R2018a kuramotoDelayed_MEX.c

load('HCP_280_IQ.mat')
load('sc_fc.mat')

% parameters
N = 82; %% number of brain regions
s = 280; % number of subjects
tSec_duration = 30; % numSeconds
stepMS = .1; % Cabral 2011: deltaT is 0.1ms
osc_freq = 40; % baseline oscillation freq.
v = 13.42; % velocity (m per sec)
w = 21.75; % global connection strength weighting


rel_localmeta40 = zeros(s,N);


%% run Kuramoto model
rng(1);
mex -R2018a kuramotoDelayed_weighted_MEX.c

% set to false to use MATLAB code instead (for checking against MEX output)
useMEX=true;


% adj matrix
for i = 1 : s
    % log of weighted connections (# of streamlines)
    net = squeeze(netsAll(:,:,i));
    a_net = mean(nonzeros(net));
    n_net = net ./ a_net;
    net = log(n_net+1);
    
    % wiring length (mm)/ mean streamlines between ROIs
    wl = squeeze(delaysAll(:,:,i));

    % execute(pure mex version)
    [theta, t, instSync, eMet, IP] = ...
                      runKuramotoDelayed_weighted(...
                      net, wl, tSec_duration, stepMS, w, v, osc_freq, useMEX);
     
    
    %%local
    local_meta40 = std(IP,0,1);
    local_meansync40 = mean(IP,1);
    local_meanSync40 = abs(local_meansync40);
    rel_localmeta40(i,:) = (local_meta40 ./ local_meanSync40)';              
                    
    
end

[r40,p40]=corr(IQ,rel_meta40,'type','Spearman');
r40_11=r40([49 12 10 3 5 14 30 7 1 43 41],:);
p40_11=p40([49 12 10 3 5 14 30 7 1 43 41],:);

[rel_r40,rel_p40]=corr(IQ,rel_localmeta40,'type','Spearman');
[max_rel_r40,index40]=max(rel_r40,[],2);
max_rel_r40_11=max_rel_r40([49 12 10 3 5 14 30 7 1 43 41],:);
index40_11=index40([49 12 10 3 5 14 30 7 1 43 41],:);


save rel_localmeta40.mat;
save rel_r40.mat;
save rel_p40.mat;
save index40_11.mat;
save max_rel_r40_11.mat;








