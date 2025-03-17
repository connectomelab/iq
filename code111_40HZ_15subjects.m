
%% Test Kuramoto (weighted) 

%% set up parameters
clear;

load('HCP_280_IQ.mat')
load('sc_fc.mat')

% parameters
s=15; %number of subjects
tSec_duration = 30; % numSeconds
stepMS = .1; % Cabral 2011: deltaT is 0.1ms
osc_freq = 40; % baseline oscillation freq.


% k loop
kappaN = 30;
kappaEnd = 30;
kappa = linspace(0.1,kappaEnd,kappaN);

velocityN = 30;
velocityEnd = 30;
velocity = linspace(0.1,velocityEnd,velocityN);

meanSync4015 = zeros(s,kappaN,velocityN);
metaStab4015 = zeros(s,kappaN,velocityN);
%Corr4015 = zeros(s,kappaN,velocityN);
MSE4015 = zeros(s,kappaN,velocityN);


%% run Kuramoto model
rng(1);
mex -R2018a kuramotoDelayed_weighted_MEX.c

% set to false to use MATLAB code instead (for checking against MEX output)
useMEX=true;


for i = 1:s
    % log of weighted connections (# of streamlines)
    net = squeeze(netsAll(:,:,i));
    a_net = mean(nonzeros(net));
    n_net = net ./ a_net;
    net = log(n_net+1);

    % empirical FC
    FC_emp = squeeze(FC(:,:,i));

    % mean streamlines between ROIs
    wl = squeeze(delaysAll(:,:,i));
    
    % k loop
    for ikappa = 1: size(kappa, 2)
             w = kappa(ikappa);
             for ivelocity = 1: size(velocity, 2)
                         v = velocity(ivelocity);


                           % execute
                          [theta, t, instSync, eMet] = ...
                                           runKuramotoDelayed_weighted(...
                                           net, wl, tSec_duration, stepMS, w, v, osc_freq, useMEX);
                          
                          meanSync4015(i,ikappa, ivelocity) = mean(instSync);
                          metaStab4015(i,ikappa, ivelocity) = eMet;


                          % simulation correlation
                          FC_sim = corr(theta); 

                          % correlation
                          %Corr4015(i,ikappa, ivelocity) = corr(FC_sim(:),FC_emp(:),'type','Spearman');
                          
                          % mean square error
                          MSE4015(i,ikappa, ivelocity) = mean((FC_sim(:)-FC_emp(:)).^2);
 
             end
             
    end
    
end

save meanSync4015.mat;
save metaStab4015.mat;
%save Corr4015.mat;
save MSE4015.mat;


MeanSync4015 = reshape(meanSync4015,[s,kappaN*velocityN]);
MetaStab4015 = reshape(metaStab4015,[s,kappaN*velocityN]);
%CORR4015 = reshape(Corr4015,[s,kappaN*velocityN]);
MSE4015re = reshape(MSE4015,[s,kappaN*velocityN]);


SmeanSync4015 = mean(MeanSync4015,1);
SmetaStab4015 = mean(MetaStab4015,1);
%SCorr4015 = mean(CORR4015,1);
MSE4015me = mean(MSE4015re,1);


R_SmeanSync4015 = reshape(SmeanSync4015,[kappaN,velocityN]);
R_SmetaStab4015 = reshape(SmetaStab4015,[kappaN,velocityN]);
%R_SCorr4015 = reshape(SCorr4015,[kappaN,velocityN]);
MSE4015end = reshape(MSE4015me,[kappaN,velocityN]);


save R_SmeanSync4015.mat;
save R_SmetaStab4015.mat;
%save R_SCorr4015.mat;
save MSE4015end.mat;





