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