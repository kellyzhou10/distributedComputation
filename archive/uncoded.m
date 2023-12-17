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
p = m;
comms = p;

Yvec = [];
for i = 1:trials
    workerTimes = 0;

    % calculate time for each worker to complete their f, take max time
    for worker = 1:p
        if (rand < pG)
            mu = muG;
        else
            mu = muB;
        end
        workerTimes = max(workerTimes, exprnd(1/mu) + tau * m / p);
    end
    Yvec = [Yvec,workerTimes];
end

disp(comms/p)

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('uncoded.mat','Yccdf','x');