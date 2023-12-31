clc
clear all
close all
format default

m = 30;
tau = 0.005;
trials = 100;
mu = 0.2;
dVec = [1,2,3,5,6,10,29]; % best not to choose d = m
SR_smallm_data = {};
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
        [time,~] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType);
        Yvec(j) = time;
    end
    [Ycdf,x] = cdfcalc(Yvec);
    Yccdf = 1 - Ycdf(1:end-1);
    semilogx(x,Yccdf);
    xlim([0.1,100]);
    if (plotType == 2)
        return;
    else
        SR_smallm_data = {"paper",x,Yccdf,-1,-1,-1,-1};
    end
end
for i = 1:length(dVec)
    d = dVec(i);
    lengths = [];
    if (plotType == 1)
        Yvec = zeros(1,trials);
        Cvec = zeros(1,trials);
        parfor j = 1:trials
            [time,~,~,comms] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
            Cvec(j) = comms;
        end
        YDvec = mean(Yvec);
        commDVec = mean(Cvec);
        addData = {string(d),commDVec,YDvec};
        SR_smallm_data = [SR_smallm_data;addData];
    elseif (plotType == 2)
        frac = 0:1:m;
        frac = frac / m;
        ACT = zeros(1,m+1);
        parfor j = 1:trials
            [~,CT,~,~] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType);
            ACT = ACT + CT;
        end
        ACT = ACT / trials;
        addData = {string(d),frac,ACT};
        SR_smallm_data = [SR_smallm_data;addData];
    elseif (plotType == 3)
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [~,~,FC,~] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType);
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
        SR_smallm_data = [SR_smallm_data;addData];
    elseif (plotType == 5)
        Yvec = zeros(1,trials);
        Cvec = zeros(1,trials);
        frac = 0:1:m;
        frac = frac / m;
        ACT = zeros(1,m+1);
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [time,CT,FC,comms] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType);
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
        SR_smallm_data = [SR_smallm_data;addData];
    end
end
if (plotType == 5)
    save('SR_smallm.mat','SR_smallm_data');
else
    graph(plotType,SR_smallm_data);
end


function [time,CT,FC,comms] = SR_smallm_subfunc(m,tau,mu,d,delta,plotType)

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
    
    j = 0;
    max_f = 0;
    prev_max_f = 0;
    CT = [];
    pool = cell(m,1);
    for k = 1:m
        pool{k} = 1:m;
    end
    while (max_f ~= m)
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
        encode = [encode;symbol];
        [max_f,~] = Intermediate_Performance(encode);
        time = workerTimes(j,2);
        if ((plotType == 2 || plotType == 3 || plotType == 5) && max_f ~= prev_max_f)
            CT(max_f) = time;
            CT(prev_max_f+1 : max_f) = time;
            prev_max_f = max_f;
        end
    end
    FC = [];
    if (plotType == 3 || plotType == 5)
        start = 1;
        for stop = 2:m
            if (CT(stop) > CT(stop-1))
                FC(floor(CT(start)*delta)+1 : floor(CT(stop)*delta)) = stop-1;
                start = stop;
            end
        end
        FC = [FC,stop];
        FC = FC / m;
    end
    CT = [0,CT];
    comms = comms / m;
end

function [max_f,exitflag] = Intermediate_Performance(encode)

    % This code is efficient only for small m (m<50).
    
    % If max_f == no_cols, the decoding is successful because this means that
    % all-one vector is in the linear span of rows of the matrix 'encode'.
    
    % As a result, this code can also be used for checking the existence of the 
    % all-one vector in the linear span of rows of the matrix 'encode'; however,
    % this code is not efficient for large m and should not be used for m>=50.
    
    options = optimoptions('intlinprog','Display','off',...
        'ObjectiveImprovementThreshold',1e-4,...
        'CutMaxIterations',25,...
        'CutGeneration','advanced',...
        'RelativeGapTolerance',1e-4);
    
    no_rows = size(encode,1);
    no_cols = size(encode,2);
    
    % no_vars = no_cols+no_rows
    
    Aeq = [eye(no_cols),-encode'];
    beq = zeros(no_cols,1);
    
    f = [-ones(1,no_cols),zeros(1,no_rows)];
    
    LB = [zeros(1,no_cols),-inf*ones(1,no_rows)];
    UB = [ones(1,no_cols),inf*ones(1,no_rows)];
    
    [~,fval,exitflag] = intlinprog(f,1:no_cols,[],[],Aeq,beq,LB,UB,[],options);
    
    if(exitflag==1)
        max_f = round(-fval);
    else
        new_encode = encode;
        new_encode(end,:) = [];
        [max_f,exitflag] = Intermediate_Performance(new_encode);
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