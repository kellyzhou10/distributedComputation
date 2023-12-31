clc
clear all
close all
format default

m = 1000;
dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,1000];
% m = 30;
% dVec = [1,2,3,5,6,10,30];
tau = 0.005;
trials = 100;
mu = 0.2;
FR_RR_data = {};
delta = 1000;

prompt = "Enter the number corresponding to the graph:" + newline + ...
    "1. communications vs time" + newline + ...
    "2. fraction completed vs average completion time" + newline + ...
    "3. time vs average function completed" + newline + ...
    "4. t vs probability of time>t (graph from research paper)" + newline + ...
    "5. save data" + newline;
plotType = input(prompt);
if (plotType == 4 || plotType == 5)
    d = 10;
    Yvec = zeros(1,trials);
    parfor j = 1:trials
        [time,~,~,comms] = FR_RR_subfunc(m,tau,mu,d,delta,plotType);
        Yvec(j) = time;
        Cvec(j) = comms;
    end
    [Ycdf,x] = cdfcalc(Yvec);
    Yccdf = 1 - Ycdf(1:end-1);
    semilogx(x,Yccdf);
    xlim([0.1,100]);
    xlabel('t');
    ylabel('Pr(T>t)');
    figure;
    [Ccdf,y] = cdfcalc(Cvec);
    Cccdf = 1 - Ccdf(1:end-1);
    plot(y,Cccdf);
    xlabel('c');
    ylabel('Pr(C>c)');
    if (plotType == 4)
        return;
    else
        FR_RR_data = {"paper",x,Yccdf,-1,-1,-1,-1};
    end
end
for i = 1:length(dVec)
    d = dVec(i);
    lengths = [];
    if (plotType == 1)
        Yvec = zeros(1,trials);
        Cvec = zeros(1,trials);
        parfor j = 1:trials
            [time,~,~,comms] = FR_RR_subfunc(m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
            Cvec(j) = comms;
        end
        YDvec = mean(Yvec);
        commDVec = mean(Cvec);
        addData = {string(d),commDVec,YDvec};
        FR_RR_data = [FR_RR_data;addData];
    elseif (plotType == 2)
        frac = 0:m/d;
        frac = frac / (m/d);
        ACT = zeros(1,m/d+1);
        parfor j = 1:trials
            [~,CT,~,~] = FR_RR_subfunc(m,tau,mu,d,delta,plotType);
            ACT = ACT + CT;
        end
        ACT = ACT / trials;
        addData = {string(d),frac,ACT};
        FR_RR_data = [FR_RR_data;addData];
    elseif (plotType == 3)
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [~,~,FC,~] = FR_RR_subfunc(m,tau,mu,d,delta,plotType);
            tempFC{j} = FC;
            lengths = [lengths,length(FC)];
        end
        len = max(lengths);
        AFC = zeros(1,len);
        for j = 1:trials
            for k = 1:len
                if (k > lengths(j))
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
        FR_RR_data = [FR_RR_data;addData];
    elseif (plotType == 5)
        Yvec = zeros(1,trials);
        Cvec = zeros(1,trials);
        frac = 0:m/d;
        frac = frac / (m/d);
        ACT = zeros(1,m/d+1);
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [time,CT,FC,comms] = FR_RR_subfunc(m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
            Cvec(j) = comms;
            ACT = ACT + CT;
            tempFC{j} = FC;
            lengths = [lengths,length(FC)];
        end
        YDvec = mean(Yvec);
        commDVec = mean(Cvec);
        ACT = ACT / trials;
        len = max(lengths);
        AFC = zeros(1,len);
        for j = 1:trials
            for k = 1:len
                if (k > lengths(j))
                    AFC(k) = AFC(k) + 1;
                else
                    AFC(k) = AFC(k) + tempFC{j}(k);
                end
            end
        end
        Tvec = 0:len-1;
        Tvec = Tvec / delta;
        AFC = AFC / trials;
        addData = {string(d),commDVec,YDvec,frac,ACT,Tvec,AFC};
        FR_RR_data = [FR_RR_data;addData];
    end
end
if (plotType == 5)
    % save('FR_RR_largem.mat');
else
    graph(plotType,FR_RR_data);
end


function [time,CT,FC,comms] = FR_RR_subfunc(m,tau,mu,d,delta,plotType)

    % define number of groups and initialize arrays: 
    %   groupTimes - holds minimum computation time for each group
    %   initialTimes - holds initial random delay value for each worker
    groups = m/d;
    groupTimes = ones(1,groups) * inf;
    initialTimes = zeros(1,m);

    for i = 1:groups

        % for each worker assigned to compute that group, calculate initial
        % random delay value and record in initialTimes; initially set the
        % minimum computation time groupT for this group to the shortest 
        % initial delay time
        groupT = inf;
        for worker = 1:d
            wTime = exprnd(1/mu);
            groupT = min(groupT,wTime);
            initialTimes(i*groups+d) = wTime;
        end
    
        % each group of workers must then compute each group of functions 
        % in a round-robin fashion
        for j = 1:groups
            % choose group number and use as index
            group = mod(i+j,groups);
            if group == 0
                group = groups;
            end
            % update groupT, the time by which the fastest worker of this 
            % group can finish computing this group of functions
            groupT = groupT + tau * d;
            % update the minimum time by which this group of functions can
            % be computed across all groups of workers
            groupTimes(group) = min(groupTimes(group),groupT);
        end
    end
    CT = [];
    FC = [];
    if (plotType == 2 || plotType == 3 || plotType == 5)
        CT = [0,sort(groupTimes)];
        if (plotType == 3 || plotType == 5)
            for i = 1:length(CT)-1
                FC(floor(CT(i)*delta)+1 : floor(CT(i+1)*delta)) = (i-1)*d;
            end
            FC = [FC,i*d];
            FC = FC / m;
        end
    end
    time = max(groupTimes);
    comms = 0;
    if (plotType == 1 || plotType == 5)
        for i = 1:m
            comms = comms + floor((time - initialTimes(i)) / tau / d);
        end
        comms = comms / m;
    end
end

function graph(input,data)
    if (input == 1)
        plot(cell2mat(data(:,2)),cell2mat(data(:,3)));
        xlabel('avg comms/worker');
        ylabel('time to complete m functions');
    else
        hold on;
        for i = 1:size(data)
            stairs(data{i,2},data{i,3});
        end
        legend(data{:,1});
        if (input == 2)
            xlabel('fraction completed');
            ylabel('average completion time');
        else
            xlabel('time');
            ylabel('average fraction completed');
        end
    end
end