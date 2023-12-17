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
dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,1000];
p = m;

YDvec = [];
commDVec = [];
for i = 1:length(dVec)
    d = dVec(i);
    Yvec = [];
    commVec = [];
    for j = 1:trials
        time = BCC_subfunc(m,tau,muG,muB,pG,d,p);
        Yvec = [Yvec,time];
        commVec = [commVec,m/d/p];
    end
    YDvec = [YDvec,mean(Yvec)];
    commDVec = [commDVec,mean(commVec)];
end

plot(commDVec,YDvec)
save('fastBCC_comms.mat','YDvec','commDVec');

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('BCC.mat','Yccdf','x');

function time = BCC_subfunc(m,tau,muG,muB,pG,d,p)

    % each worker completes all batches in random order
    % take the fastest completion time for each batch across workers
    % take the max time of all batches
    
    batches = m/d;
    batchTimes = ones(batches,1) * inf;
    for worker = 1:p
        if (rand < pG)
            mu = muG;
        else
            mu = muB;
        end
        workerTime = exprnd(1/mu);
    
        % for each randomly selected batch in the remaining batches,
        % compute the total time at which this worker will finish the
        % batch, and record the minimum of this time and the previously
        % recorded batchTime of this batch
    
        I = 1:batches;
        for batch = 1:batches
            coupon = randsrc(1,1,I);
            I = setdiff(I,coupon);
            workerTime = workerTime + tau * d;
            batchTimes(coupon) = min(batchTimes(coupon), workerTime);
        end
    end
    
    time = max(batchTimes);

end