clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 100;
muG = 0.4;
muB = 0.2;
pG = 0;
d = 10;
p = m;
comms = m/d;

Yvec = [];
for i = 1:trials
    
    % take the max time out of all groups
    groupTimes = zeros(1,m/d);
    for group = 1:m/d

        % take the time for the fastest worker in each group
        Ysub = 1000;
        for worker = 1:p/m*d
            if (rand < pG)
                mu = muG;
            else
                mu = muB;
            end
            Ysub = min(Ysub, exprnd(1/mu) + tau * d);
        end
        groupTimes(group) = Ysub;
    end
    Yvec = [Yvec,max(groupTimes)];
end

disp(comms/p)

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('FR.mat','Yccdf','x');