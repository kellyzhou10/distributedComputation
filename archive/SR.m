clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 1;
muG = 0.4;
muB = 0.2;
pG = 0;
d = 10;
k = m;
p = m;
global comms;
comms = 0;

Yvec = [];
for i = 1:trials
    time = SR_subfunc(m,tau,muG,muB,pG,d,k,p);
    Yvec = [Yvec,time];
end

disp(comms/p)

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('SR.mat','Yccdf','x');

function time = SR_subfunc(m,tau,muG,muB,pG,d,k,p)

    % for each worker, generate the times to compute each symbol, then
    % generate a row for each symbol using the fastest workers, adding to
    % the matrix until it is solvable
    
    global comms;
    encode = [];
    workerTimes = zeros(p,ceil(1+(p-1)/d));
    
    % setting time of computation for each worker; first computation is a
    % source symbol
    for worker = 1:p
        if (rand < pG)
            mu = muG;
        else
            mu = muB;
        end
        workerTimes(worker,1) = exprnd(1/mu) + tau * m / p;
        for j = 2:ceil(1+(p-1)/d)
            workerTimes(worker,j) = tau * m / p * d + workerTimes(worker,j-1);
        end
    end
    
    % while matrix is unsolvable using Gaussian elim, generate encoded
    % symbols corresponding to the generated times
    j = 0;
    temp = [ones(p,1),workerTimes(:,1)];
    while (true)
        [time,worker] = min(temp(:,2));
    
        % create new row for every symbol; each worker's first is a source
        % symbol, all others are encoded
        symbol = zeros(1,p);
        if (temp(worker,1) == 1)
            symbol(worker) = 1;
        else
            I = 1:p;
            I = setdiff(I,worker);
            for counter = 1:d
                f = randsrc(1,1,I);
                I = setdiff(I,f);
                symbol(f) = 1;
            end
        end
        comms = comms + 1;
        encode = [encode;symbol];
        if (j >= p && rank(encode) == p)
            break;
        end
        temp(worker,1) = temp(worker,1) + 1;
        temp(worker,2) = workerTimes(worker,temp(worker,1));
        j = j + 1;
    end
end