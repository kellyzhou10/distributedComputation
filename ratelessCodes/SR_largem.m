clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 1;
mu = 0.2;
dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,999]; % best not to choose d = m
SR_largem_data = {};
delta = 1000;

prompt = "Enter the number corresponding to the graph:" + newline + ...
    "1. communications vs time" + newline + ...
    "2. t vs probability of time>t (graph from research paper)" + newline + ...
    "3. save data" + newline;
plotType = input(prompt);
if (plotType == 2 || plotType == 3)
    d = 10;
    Yvec = zeros(1,trials);
    parfor j = 1:trials
        [time,~] = SR_largem_subfunc(m,tau,mu,d);
        Yvec(j) = time;
    end
    [Ycdf,x] = cdfcalc(Yvec);
    Yccdf = 1 - Ycdf(1:end-1);
    semilogx(x,Yccdf);
    xlim([0.1,100]);
    if (plotType == 2)
        return;
    else
        SR_largem_data = {"paper",x,Yccdf,-1,-1,-1,-1};
    end
end
for i = 1:length(dVec)
    d = dVec(i);
    Yvec = zeros(1,trials);
    Cvec = zeros(1,trials);
    parfor j = 1:trials
        [time,comms] = SR_largem_subfunc(m,tau,mu,d);
        Yvec(j) = time;
        Cvec(j) = comms;
    end
    YDvec = mean(Yvec);
    commDVec = mean(Cvec);
    addData = {string(d),commDVec,YDvec,-1,-1,-1,-1};
    SR_largem_data = [SR_largem_data;addData];
end
if (plotType == 3)
    save('SR_largem.mat','SR_largem_data');
else
    graph(plotType,SR_largem_data);
end


function [time,comms] = SR_largem_subfunc(m,tau,mu,d)

    % for each worker, generate the times to compute each symbol, then
    % generate a row for each symbol using the fastest workers, adding to
    % the matrix until it is solvable
    
    comms = 0;
    encode = [];
    workerTimes = zeros(m,3);
    
    % setting time of computation for each worker; first computation is a
    % source symbol
    for worker = 1:m
        workerTimes(worker,1) = worker;
        workerTimes(worker,2) = exprnd(1/mu) + tau;
        workerTimes(worker,3) = 1;
    end
    lastMin = min(workerTimes(:,2));
    firstMax = max(workerTimes(:,2));
    round = 0;
    while (lastMin < firstMax)
        tempTimes = zeros(m,3);
        for worker = 1:m
            tempTimes(worker,1) = worker;
            tempTimes(worker,2) = tau * d + workerTimes(round*m+worker,2);
        end
        round = round + 1;
        lastMin = min(tempTimes(:,2));
        workerTimes = [workerTimes;tempTimes];
    end
    workerTimes = sortrows(workerTimes,2);
    
    % while matrix is unsolvable, generate encoded symbols corresponding to 
    % the generated times
    j = 0;
    prevRank = 0;
    onesVec = ones(1,m);
    allOnes = false;
    pool = cell(m,1);
    for k = 1:m
        pool{k} = 1:m;
    end
    while (~allOnes)
        j = j + 1;

        % create new row for every symbol; each worker's first is a source
        % symbol, all others are encoded
        symbol = zeros(1,m);
        worker = workerTimes(j,1);
        if (workerTimes(j,3) == 1)
            symbol(worker) = 1;
            pool{worker} = setdiff(pool{worker},worker);
        else
            if (length(pool{worker}) < d)
                continue;
            end
            for counter = 1:d
                f = randsrc(1,1,pool{worker});
                pool{worker} = setdiff(pool{worker},f);
                symbol(f) = 1;
            end
        end
        comms = comms + 1;
        curRank = rank([encode;symbol]);
        if (curRank > prevRank)
            prevRank = curRank;
            encode = [encode;symbol];
            if (rank([encode;onesVec]) == curRank)
                allOnes = true;
            end
        end
    end
    time = workerTimes(j,2);
    comms = comms / m;
end

function graph(input,data)
    if (input == 1)
        plot(cell2mat(data(:,2)),cell2mat(data(:,3)));
    else
        hold on;
        for i = 1:size(data)
            stairs(data{i,2},data{i,3});
        end
        legend(data{:,1});
    end
end