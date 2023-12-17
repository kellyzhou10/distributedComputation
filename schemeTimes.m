clc
clear all
close all

load('allAtOne.mat');
semilogx(x,Yccdf);

xlim([0.1,100]);
xlabel('t');
ylabel('Pr(T>t)');
hold on

load('uncoded.mat');
semilogx(x,Yccdf);

load('fastFR_RR.mat');
semilogx(x,Yccdf);

load('BCC.mat');
semilogx(x,Yccdf);

load('fastLT.mat');
semilogx(x,Yccdf);

load('SR.mat');
semilogx(x,Yccdf);

legend('One','Unc','Rep','BCC','LT','SR');