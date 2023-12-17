clc
clear all
close all
format default

m = 1000;
me = 2000;
tau = 0.005;
trials = 1;
muG = 0.4;
muB = 0.2;
pG = 0;
p = m;
global comms;
comms = 0;
% c = 0.36
% K = m
% R = 79
% delta = 0.95

Yvec = [];
for i = 1:trials
    
    % for each worker, pick a random number d according to robust 
    % distribution of random fs & compute the total time to complete them;
    % repeat immediately for a total of me/m rounds;
    % while the rank of encode is not m, generate encoded symbols
    % corresponding with the fastest times

    global comms;
    encode = [];
    workerTimes = zeros(p*me/m,2);
    for worker = 1:p
        if (rand < pG)
            mu = muG;
        else
            mu = muB;
        end
        workerTime = exprnd(1/mu);
        for round = 0:me/m-1
            % setting degree and time of computation
            d = randsrc(1,1,[1 2 3 4 5 6 7 8 9 10 11 12; ...
                0.0790 0.0395 0.0263 0.0198 0.0158 0.0132 ...
                0.0113 0.0099 0.0088 0.0079 0.0072 0.7613]);
            workerTimes(worker+p*round,1) = d;
            workerTimes(worker+p*round,2) = workerTime + tau * d * m / p;
        end
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
        comms = comms + 1;
        encode = [encode;symbol];
    end
    Yvec = [Yvec,workerTimes(j,2)];
end

disp(comms/p)

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
semilogx(x,Yccdf);
xlim([0.1,100]);