clc
clear all
close all

load('FR_RR_largem.mat');
load('FR_rand_largem.mat');
load('BCC_largem.mat');
load('LT_largem.mat');
load('SR_largem.mat');

prompt1 = "Enter the number corresponding to the graph:" + newline + ...
    "1. communications vs time" + newline + ...
    "2. fraction completed vs average completion time" + newline + ...
    "3. time vs average function completed" + newline + ...
    "4. t vs probability of time>t (graph from research paper)" + newline;
graphType = input(prompt1);
if (graphType == 1)
    I = 2:16;
    plot(cell2mat(FR_RR_data(I,2)),cell2mat(FR_RR_data(I,3)));
    hold on;
    plot(cell2mat(FR_rand_data(I,2)),cell2mat(FR_rand_data(I,3)));
    plot(cell2mat(BCC_data(I,2)),cell2mat(BCC_data(I,3)));
    plot(cell2mat(LT_largem_data(I,2)),cell2mat(LT_largem_data(I,3)));
    plot(cell2mat(SR_largem_data(I,2)),cell2mat(SR_largem_data(I,3)));
    xlabel('avg comms/worker');
    ylabel('time to complete m functions');
elseif (graphType == 4)
    semilogx(FR_RR_data{1,2},FR_RR_data{1,3});
    hold on;
    semilogx(FR_rand_data{1,2},FR_rand_data{1,3});
    semilogx(BCC_data{1,2},BCC_data{1,3});
    semilogx(LT_largem_data{1,2},LT_largem_data{1,3});
    semilogx(SR_largem_data{1,2},SR_largem_data{1,3});
    xlim([0.1,100]);
    xlabel('t');
    ylabel('Pr(T>t)');
else
    load('FR_RR_smallm.mat');
    load('FR_rand_smallm.mat');
    load('BCC_smallm.mat');
    load('LT_smallm.mat');
    load('SR_smallm.mat');
    prompt2 = "Select d in: [1,2,3,5,6,10,30];" + newline;
    d = input(prompt2);
    for i = 1:length(FR_rand_data)
        if (str2num(FR_rand_data{i,1}) == d)
            break;
        end
    end
    if (graphType == 2)
        x = 4;
        y = 5;
        xlabel('fraction completed');
        ylabel('average completion time');
    else
        x = 6;
        y = 7;
        xlabel('time');
        ylabel('average fraction completed');
    end
    hold on;
    stairs(FR_RR_data{i,x},FR_RR_data{i,y});
    stairs(FR_rand_data{i,x},FR_rand_data{i,y});
    stairs(BCC_data{i,x},BCC_data{i,y});
    stairs(LT_smallm_data{i,x},LT_smallm_data{i,y});
    stairs(SR_smallm_data{i,x},SR_smallm_data{i,y});
end
legend('FR RR','FR rand','BCC','LT','SR');