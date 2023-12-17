clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 100;
mu = 0.2;
dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,1000];
data = {};
delta = 1000;

prompt = "Enter the number corresponding to the graph:" + newline + ...
    "1. communications vs time" + newline + ...
    "2. fraction completed vs average completion time" + newline + ...
    "3. time vs average function completed" + newline + ...
    "4. t vs probability of time>t (graph from research paper)" + newline;
plotType = input(prompt);
for i = 1:length(dVec)
    d = dVec(i);
    lengths = [];
    if (plotType == 1)
        Yvec = zeros(1,trials);
        Cvec = zeros(1,trials);
        parfor j = 1:trials
            [time,~,~,comms] = FR_subfunc(m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
            Cvec(j) = comms;
        end
        YDvec = mean(Yvec);
        commDVec = mean(Cvec);
        addData = {string(d),commDVec,YDvec};
        data = [data;addData];
    elseif (plotType == 2)
        frac = 0:m/d;
        frac = frac / (m/d);
        ACT = zeros(1,m/d+1);
        parfor j = 1:trials
            [~,CT,~,~] = FR_subfunc(m,tau,mu,d,delta,plotType);
            ACT = ACT + CT;
        end
        ACT = ACT / trials;
        addData = {string(d),frac,ACT};
        data = [data;addData];
    elseif (plotType == 3)
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [~,~,FC,~] = FR_subfunc(m,tau,mu,d,delta,plotType);
            tempFC{j} = FC;
            lengths = [lengths,length(FC)];
        end
        len = max(lengths);
        AFC = zeros(1,len);
        for j = 1:trials
            for k = 1:len
                if k > lengths(j)
                    AFC(k) = AFC(k) + 1;
                else
                    AFC(k) = AFC(k) + tempFC{j}(k);
                end
            end
        end
        Tvec = 0:len-1;
        Tvec = Tvec / delta;
        AFC = AFC / trials;
        addData = {string(d),Tvec,AFC};
        data = [data;addData];
    else
        d = 10;
        Yvec = zeros(1,trials);
        parfor j = 1:trials
            [time,~,~,~] = FR_subfunc(m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
        end
        [Ycdf,x] = cdfcalc(Yvec);
        Yccdf = 1 - Ycdf(1:end-1);
        semilogx(x,Yccdf);
        xlim([0.1,100]);
        return;
    end
end
graph(plotType,data);


function [time,CT,FC,comms] = FR_subfunc(m,tau,mu,d,delta,plotType)

    % each group of workers computes their own group of functions and
    % then compute the next group, repeat until all groups are covered;
    % computation is complete when all groups are covered by the worker
    % with the shortest time

    groups = m/d;
    groupTimes = ones(1,groups) * inf;
    initialTimes = zeros(1,m);
    for i = 1:groups
        groupT = inf;
        for worker = 1:d
            wTime = exprnd(1/mu);
            groupT = min(groupT,wTime + tau * d);
            initialTimes(i*groups+d) = wTime;
        end
        groupTimes(i) = groupT;
    end
    CT = [];
    FC = [];
    if plotType == 2 || plotType == 3
        CT = [0,sort(groupTimes)];
        if plotType == 3
            for i = 1:length(CT)-1
                FC(floor(CT(i)*delta)+1 : floor(CT(i+1)*delta)) = (i-1)*d;
            end
            FC = [FC,i*d];
            FC = FC / m;
        end
    end
    time = max(groupTimes);
    comms = 0;
    for i = 1:m
        comms = comms + floor((time - initialTimes(i)) / tau / d);
    end
    comms = comms / m;
end

function graph(input,data)
    if input == 1
        plot(cell2mat(data(:,2)),cell2mat(data(:,3)));
    else
        hold on;
        for i = 1:size(data)
            stairs(data{i,2},data{i,3});
        end
        legend(data{:,1});
    end
end