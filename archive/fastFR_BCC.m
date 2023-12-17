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
        time = fastFR_BCC_subfunc(m,tau,muG,muB,pG,d,p);
        Yvec = [Yvec,time];
        commVec = [commVec,m/d/p];
    end
    YDvec = [YDvec,mean(Yvec)];
    commDVec = [commDVec,mean(commVec)];
end

plot(commDVec,YDvec)
save('fastFR_BCC_comms.mat','YDvec','commDVec');

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
save('fastFR_BCC.mat','Yccdf','x');


function time = fastFR_BCC_subfunc(m,tau,muG,muB,pG,d,p)

    % each group of workers computes their own group of functions and
    % then compute other groups, chosen at random; computation is complete
    % when all groups are covered by the worker with the shortest time
    
    groups = m/d;
    groupTimes = ones(groups,1) * inf;
    for group = 1:groups
        groupTime = inf;
        for worker = 1:floor(p/groups)
            if (rand < pG)
                mu = muG;
            else
                mu = muB;
            end
            groupTime = min(groupTime,exprnd(1/mu));
        end
        groupTime = groupTime + tau * d;
        groupTimes(group) = min(groupTimes(group),groupTime);
    
        % for each randomly selected group in the remaining groups,
        % compute the total time at which this worker will finish the
        % group, and record the minimum of this time and the previously
        % recorded groupTime of this group
    
        I = 1:groups;
        I = setdiff(I,group);
        for remGroup = 1:groups-1
            coupon = randsrc(1,1,I);
            I = setdiff(I,coupon);
            groupTime = groupTime + tau * d;
            groupTimes(coupon) = min(groupTimes(coupon), groupTime);
        end
    end
    
    time = max(groupTimes);

end