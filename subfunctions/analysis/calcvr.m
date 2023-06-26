function distMat = calcvr(spTimes, tau) %, nTargets, tempnum)
% distMat = calcdistmat(spTimes, tau)
% Both spTimes and tau are in seconds, not ms (as calc_distance uses)
% This function runs faster if the mex file calcvrHelper.c is compiled.

% this barfs if there is a trial with no spikes

if exist('calcvrHelper', 'file') == 3
	yFun = @calcvrHelper;
else
	yFun = @calcvrHelperSlow;
end
if exist('calcvrSorter', 'file') == 3
	ySor = @calcvrSorter;
else
	ySor = @calcvrSorterSlow;
end
if nargin == 3
	error('Number of inputs must be 2 or 4.')
end

numTaus = length(tau); %get # of taus to be used
numTrains = numel(spTimes); %spTimes is a cell containing spike trains in each cell
distMat = zeros(numTrains, numTrains, numTaus); 
% distance matrix: the VR of [numTrains] trains at different tau=[numTaus]
mylens = cellfun(@length,spTimes); 
% applies length() to each element (spike train times) in spTimes
for ind1 = 1:numTrains
	spTimesPart1 = spTimes{ind1};
	for ind2 = ind1 + 1:numTrains
		spTimesPart2 = spTimes{ind2};
		if mylens(ind1) + mylens(ind2) ~= 0 
            [t,w] = ySor(spTimesPart1,spTimesPart2); % calculate difference function
			minusDifft = -diff(t); % convert time stamps to negative intervals
			for tauNum = 1:numTaus
                % create decaying exponential kernel
				expMinusDifft = exp(minusDifft/tau(tauNum));%=exp(-t/tau)
                
                % calculate VR distance using integral of w convolved with
                % decaying exponential function
				distMat(ind1,ind2,tauNum) = ...
                    sum((yFun(w, expMinusDifft) + w) .^2 .* ([1 - expMinusDifft.*expMinusDifft; 1]))*0.5;
				distMat(ind2,ind1,tauNum) = distMat(ind1,ind2,tauNum);
			end
		end
	end
end

function y = calcvrHelperSlow(w, e)
y = zeros(numel(w),1);
for ii = 1:numel(w) - 1
	y(ii + 1) = (y(ii) + w(ii))*e(ii);
end

function [tout,wout] = calcvrSorterSlow(t1,t2)

% tout: time of all events between both trains
% wout: the difference function between both trains

n1 = length(t1);
n2 = length(t2);
nout = n1+n2;

tout = zeros(nout, 1);
wout = zeros(nout, 1);

c1=1;
c2=1;
if(isempty(t1))
    t1 = inf;
end
if(isempty(t2))
    t2 = inf;
end
for count=1:nout
    if(t1(c1) < t2(c2))
        tout(count) = t1(c1);
        wout(count) = 1;
        c1 = c1+1;
        if(c1>n1)
            break;
        end
    else
        tout(count) = t2(c2);
        wout(count) = -1;
        c2 = c2+1;
        if(c2>n2)
            break;
        end
    end
end
count=count+1;
if(c1>n1)
    if(c2<=n2)
        tout(count:nout) = t2(c2:n2);
        wout(count:nout) = -1;
    end
elseif(c2>n2)
    if(c1<=n1)
        tout(count:nout) = t1(c1:n1);
        wout(count:nout) = 1;
    end
end

