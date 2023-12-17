clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 4;
muG = 0.4;
muB = 0.2;
pG = 0;
p = m;
% c = 0.36
% K = m
% R = 79
% delta = 0.95

Yvec = [];
for i = 1:trials
    
    % for each worker, pick a random number d according to robust 
    % distribution of random fs & compute the total time to complete them;
    % workers do not wait and will receive another set of d fs to compute
    % as soon as it is finished with one set (so we must generate enough
    % times to complete each set for each worker to determine the minimum
    % time to complete the entire summation);
    % while the rank of encode is not m, generate encoded symbols
    % corresponding with the fastest times

    encode = [];
    workerTimes = zeros(p,2);
    for worker = 1:p

        % setting degree and time of computation
        d = randsrc(1,1,[1 2 3 4 5 6 7 8 9 10 11 12; ...
            0.0790 0.0395 0.0263 0.0198 0.0158 0.0132 ...
            0.0113 0.0099 0.0088 0.0079 0.0072 0.7613]);
        workerTimes(worker,1) = d;
        if (rand < pG)
            mu = muG;
        else
            mu = muB;
        end
        workerTimes(worker,2) = exprnd(1/mu) + tau * d * m / p;
    end
    firstMax = max(workerTimes(:,2));
    lastMin = 0;

    % while the minimum of the current round of computations across all 
    % workers is greater than the maximum of the first round
    round = 0;
    while (lastMin < firstMax)

        tempTimes = zeros(p,2);
        for worker = 1:p

            % setting degree and time of computation
            d = randsrc(1,1,[1 2 3 4 5 6 7 8 9 10 11 12; ...
                0.0790 0.0395 0.0263 0.0198 0.0158 0.0132 ...
                0.0113 0.0099 0.0088 0.0079 0.0072 0.7613]);
            tempTimes(worker,1) = d;
            tempTimes(worker,2) = tau * d * m / p + workerTimes(round*p+worker,2);
        end
        round = round + 1;
        lastMin = min(tempTimes(:,2));
        workerTimes = [workerTimes;tempTimes];
    end
    workerTimes = sortrows(workerTimes,2);

    % while matrix is unsolvable using Gaussian elim, generate encoded 
    % symbols corresponding to the generated times
    j = 0;
    while (j < p || rank(encode) ~= p)
        j = j + 1;

        % create new row
        I = 1:p;
        symbol = zeros(1,p);
        for counter = 1:workerTimes(j,1)
            f = randsrc(1,1,I);
            I = setdiff(I,f);
            symbol(f) = 1;
        end
        encode = [encode;symbol];
    end
    Yvec = [Yvec,workerTimes(j,2)];
end

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('fastLT.mat','Yccdf','x');