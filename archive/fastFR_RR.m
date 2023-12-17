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
%dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,1000];
dVec = [1,2,4,5,8,10,20,40,50,100];

p = 10;
YDvec = [];
commDVec = [];
for i = 1:length(dVec)
    d = dVec(i);
    Yvec = [];
    commVec = [];
    for j = 1:trials
        time = fastFR_RR_subfunc(m,tau,muG,muB,pG,d,p);
        Yvec = [Yvec,time];
        commVec = [commVec,m/d/p];
    end
    YDvec = [YDvec,mean(Yvec)];
    commDVec = [commDVec,mean(commVec)];
end
disp(YDvec)
disp(commDVec)

hold on
plot(commDVec,YDvec)


%save('fastFR_RR_comms.mat','YDvec','commDVec');

% d = 1000;
% YDvec = [];
% commDVec = [];
% for i = 1:length(pVec)
%     p = pVec(i);
%     Yvec = [];
%     commVec = [];
%     for j = 1:trials
%         time = fastFR_RR_subfunc(m,tau,muG,muB,pG,d,p);
%         Yvec = [Yvec,time];
%         commVec = [commVec,m/d/p];
%     end
%     YDvec = [YDvec,mean(Yvec)];
%     commDVec = [commDVec,mean(commVec)];
% end
% 
% plot(commDVec,YDvec)
% save('fastFR_RR_comms.mat','YDvec','commDVec');

[Ycdf,x] = cdfcalc(Yvec);
Yccdf = 1 - Ycdf(1:end-1);
%save('fastFR_RR.mat','Yccdf','x');


function time = fastFR_RR_subfunc(m,tau,muG,muB,pG,d,p)

    % each group of workers computes their own group of functions and
    % then compute the next group, repeat until all groups are covered;
    % computation is complete when all groups are covered by the worker
    % with the shortest time

    groups = m/d;
    groupTimes = ones(groups,1) * inf;
    for j = 1:groups
        groupTime = inf;
        for worker = 1:floor(p/groups)
            if (rand < pG)
                mu = muG;
            else
                mu = muB;
            end
            groupTime = min(groupTime,exprnd(1/mu));
        end
    
        % for each remaining group, compute the total time at which this
        % worker will finish the group, and record the minimum of this time
        % and the previously recorded groupTime of this group
    
        for k = 1:groups
            group = mod(j+k-1,groups);
            if group == 0
                group = groups;
            end
            groupTime = groupTime + tau * d;
            groupTimes(group) = min(groupTimes(group),groupTime);
        end
    end
    
    time = max(groupTimes);

end