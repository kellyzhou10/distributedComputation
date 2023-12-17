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


