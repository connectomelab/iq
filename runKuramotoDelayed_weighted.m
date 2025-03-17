function [theta_sin,t,instSync,eMet] = ...
    runKuramotoDelayed_weighted(net, e_dist, tSec, step, K, V, osc, useMEX)


n = size(net,1);
p.N = n;
p.M = net;
p.K = single(K); % global weighting

% oscillation freq, per node, with noise (1%)
natOsc = single(2*pi*osc);
noise = randn(n,1) * 0.01 * natOsc;
natOsc = natOsc + noise;
p.osc = natOsc; 

nItersPerSec = single(1000 / step); % step=1 is one millisecond;

% 1 metre per sec is the same speed as 1mm per millisec.
p.D = ceil((e_dist ./ V) ./ step); % p.D = DELAY (in terms of numbers of steps) 
% delays - sets delays to zero for non-connections
p.D = p.D .* (double(p.M>0));

% initialise theta - evenly spaced-out values
thetaInit = single(((1:p.N) ./ p.N) * pi*2);

tic;
if useMEX
    % MEX VERSION:
    [theta] = kuramotoDelayed_weighted_MEX(...
        single(thetaInit(:)),...
        double(tSec), ...
        double(nItersPerSec),...
        sparse(double(p.M)),... % sparse MUST be double
        sparse(double(p.D)),... % sparse MUST be double
        single(p.K), ...
        single(p.osc),...
        single(max(p.D(:))));

    % remove buffer, and transpose (IMPORTANT)
    theta = theta(:,(max(p.D(:))+1):end)';
else
    fprintf('WARN: Using MATLAB version of kuramoto \n');
    % MATLAB VERSION:
    [theta] = kuramotoDelayed_weighted(thetaInit(:),tSec,nItersPerSec,p);
    % Faster version:
    %[theta] = kuramotoDelayed_weighted_faster(thetaInit(:),tSec,nItersPerSec,p);
end
toc;

% transpose
%theta = theta';

totTime = tSec*nItersPerSec;
% time in seconds
t = (1:totTime) ./ nItersPerSec;


%if ~exist('dontCrop', 'var') || ~dontCrop
    % remove first 10 %
    keepI = floor(.1*totTime)+1;
    t = t(keepI:end);
    theta = theta(keepI:end,:);
%end

% speed things up (downsample to ms resolution)
sparseKeep = (1:((1/step)):size(theta,1));
theta = theta(sparseKeep,:);
t = t(sparseKeep);

% sync and meta:
IP = exp(theta*1i);

% instantaneous synchrony
eCoher = (1/p.N)*sum(IP,2);
instSync = abs(eCoher);

% metastability:
eMet = std(instSync);

% return sin of theta
theta_sin = sin(theta);

end
