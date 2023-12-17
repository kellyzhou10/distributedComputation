clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 10;
muG = 0.4;
muB = 0.2;
pG = 0;
comms = 1;

Yvec = [];
for i = 1:trials

    % calculate time for 1 worker to complete all fs
    if (rand < pG)
        mu = muG;
    else
        mu = muB;
    end
    Yj = exprnd(1/mu) + tau * m;
    Yvec = [Yvec,Yj];
end

disp(comms)

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('allAtOne.mat','Yccdf','x')