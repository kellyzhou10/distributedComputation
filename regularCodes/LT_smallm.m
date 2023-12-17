clc
clear all
close all
format default

m = 30;
me = 1200;
tau = 0.005;
trials = 4;
mu = 0.2;
dVec = [1,2,3,5,6,10,30];
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
            [time,~,~,comms] = LT_subfunc(me,m,tau,mu,d,delta,plotType);
            Yvec(j) = time;
            Cvec(j) = comms;
        end
        YDvec = mean(Yvec);
        commDVec = mean(Cvec);
        addData = {string(d),commDVec,YDvec};
        data = [data;addData];
    elseif (plotType == 2)
        frac = 0:1:m;
        frac = frac / m;
        ACT = zeros(1,m+1);
        parfor j = 1:trials
            [~,CT,~,~] = LT_subfunc(me,m,tau,mu,d,delta,plotType);
            ACT = ACT + CT;
        end
        ACT = ACT / trials;
        addData = {string(d),frac,ACT};
        data = [data;addData];
    elseif (plotType == 3)
        tempFC = cell(trials,1);
        parfor j = 1:trials
            [~,~,FC,~] = LT_subfunc(me,m,tau,mu,d,delta,plotType);
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
            [time,~,~,~] = LT_subfunc(me,m,tau,mu,d,delta,plotType);
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

function [time,CT,FC,comms] = LT_subfunc(me,m,tau,mu,d,delta,plotType)

    % for each worker, pick a random number d according to robust 
    % distribution of random fs & compute the total time to complete them;
    % repeat immediately for a total of me/m rounds;
    % while the rank of encode is not m, generate encoded symbols
    % corresponding with the fastest times

    comms = 0;
    encode = [];
    workerTimes = zeros(me,2);
    [deg_vec,p_vec,~] = RobustSolitonDistribution(m,d);
    for worker = 1:m
        workerTime = exprnd(1/mu);
        for round = 0:me/m-1
            % setting degree and time of computation
            deg = randsrc(1,1,[deg_vec;p_vec]);
            workerTimes(worker+m*round,1) = deg;
            if round == 0
                workerTimes(worker,2) = workerTime + tau * deg;
            else
                workerTimes(worker+m*round,2) = workerTimes(worker+m*(round-1),2) + tau * deg;
            end
        end
    end
    workerTimes = sortrows(workerTimes,2);

    % while matrix is unsolvable using Gaussian elim, generate encoded 
    % symbols corresponding to the generated times
    j = 0;
    max_f = 0;
    prev_max_f = 0;
    CT = [];
    while (max_f ~= m)
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
        encode = [encode;symbol];
        [max_f,~] = Intermediate_Performance(encode);
        time = workerTimes(j,2);
        if (plotType == 2 || plotType == 3) && max_f ~= prev_max_f
            CT(max_f) = time;
            CT(prev_max_f+1 : max_f) = time;
            prev_max_f = max_f;
        end
    end
    FC = [];
    if plotType == 3
        start = 1;
        for stop = 2:m
            if CT(stop) > CT(stop-1)
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