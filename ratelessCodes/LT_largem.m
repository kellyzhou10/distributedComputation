clc
clear all
close all
format default

m = 1000;
tau = 0.005;
trials = 100;
mu = 0.2;
dVec = [1,2,4,5,8,10,20,40,50,100,125,200,250,500,1000];
LT_largem_data = {};
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
        [time,comms] = LT_largem_subfunc(m,tau,mu,d);
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
    if (plotType == 2)
        return;
    else
        LT_largem_data = {"paper",x,Yccdf,-1,-1,-1,-1};
    end
end
for i = 1:length(dVec)
    d = dVec(i);
    Yvec = zeros(1,trials);
    Cvec = zeros(1,trials);
    parfor j = 1:trials
        [time,comms] = LT_largem_subfunc(m,tau,mu,d);
        Yvec(j) = time;
        Cvec(j) = comms;
    end
    YDvec = mean(Yvec);
    commDVec = mean(Cvec);
    addData = {string(d),commDVec,YDvec,-1,-1,-1,-1};
    LT_largem_data = [LT_largem_data;addData];
end
if (plotType == 3)
    % save('LT_largem.mat');
else
    graph(plotType,LT_largem_data);
end


function [deg_vec,p_vec,opt_dave] = RobustSolitonDistribution(m,desired_dave)
    
    % m is the total number of files
    % desired_dave is the desired average degree d
    
    delta = 0.95;
    min_diff = inf;
    for c = 0.001:0.001:1
        R = c*log(m/delta)*sqrt(m);
        if(R>=1)
            p = zeros(1,m);
            for i = 1:floor(m/R)-1
                p(i) = R/(i*m);
            end
            p(floor(m/R)) = 1-sum(p);
            dave = sum((1:m).*p);
            if(abs(dave-desired_dave)<min_diff)
                min_diff = abs(dave-desired_dave);
                opt_p = p;
                opt_dave = dave;
            end
        end
    end
    
    deg_vec = find(opt_p~=0); 
    % deg_vec is the vector of degrees d among which each worker chooses one
    % according to the probabilities specified in p_vec
    
    p_vec = opt_p(deg_vec);
    % p_vec is the vector of probabilities associated with degrees d in deg_vec

end

function [time,comms] = LT_largem_subfunc(m,tau,mu,d)

    % for each worker, pick a random number d according to robust 
    % distribution of random fs & compute the total time to complete them;
    % workers do not wait and will receive another set of d fs to compute
    % as soon as it is finished with one set (so we must generate enough
    % times to complete each set for each worker to determine the minimum
    % time to complete the entire summation);
    % while the rank of encode is not m, generate encoded symbols
    % corresponding with the fastest times
    
    comms = 0;
    encode = [];
    workerTimes = zeros(m,2);
    [deg_vec,p_vec,~] = RobustSolitonDistribution(m,d);
    for worker = 1:m

        % setting degree and time of computation
        deg = randsrc(1,1,[deg_vec; p_vec]);
        workerTimes(worker,1) = deg;
        workerTimes(worker,2) = exprnd(1/mu) + tau * deg;
    end
    firstMax = max(workerTimes(:,2));
    lastMin = 0;

    % while the minimum of the current round of computations across all 
    % workers is greater than the maximum of the first round
    round = 0;
    while (lastMin < firstMax)

        tempTimes = zeros(m,2);
        for worker = 1:m

            % setting degree and time of computation
            deg = randsrc(1,1,[deg_vec; p_vec]);
            tempTimes(worker,1) = deg;
            tempTimes(worker,2) = tau * deg + workerTimes(round*m+worker,2);
        end
        round = round + 1;
        lastMin = min(tempTimes(:,2));
        workerTimes = [workerTimes;tempTimes];
    end
    workerTimes = sortrows(workerTimes,2);

    % while matrix is unsolvable using Gaussian elim, generate encoded 
    % symbols corresponding to the generated times
    j = 0;
    prevRank = 0;
    onesVec = ones(1,m);
    allOnes = false;
    while (~allOnes)
        j = j + 1;

        % create new row
        I = 1:m;
        symbol = zeros(1,m);
        for counter = 1:workerTimes(j,1)
            f = randsrc(1,1,I);
            I = setdiff(I,f);
            symbol(f) = 1;
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