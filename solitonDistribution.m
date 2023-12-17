% function [d,P] = solitonDistribution(K,delta,c)
    delta = 0.95;
    K = 1000;
    for c = 0.36
        R = round(c*log(K/delta)*sqrt(K))
        p = zeros(1,K);
        for i = 1:floor(K/R)-1
            p(i) = R/(i*K);
        end
        p(floor(K/R)) = 1-sum(p);
        ave = sum((1:K).*p);
        disp(p(1:12));
    end
% end