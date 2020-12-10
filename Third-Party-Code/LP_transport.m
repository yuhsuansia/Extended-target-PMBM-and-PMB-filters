function [Cmin,q] = LP_transport(C,qh,qj)

[H,N] = size(C);
c = reshape(C',H*N,1);
beq = [qh;qj];
Aeq = sparse([kron(eye(H),ones(1,N));repmat(eye(N),1,H)]);

%MATLAB solver
c(isinf(c)) = max(c(~isinf(c)))+1e6;
options = optimoptions('linprog','Display','off');
[x,Cmin] = linprog(c,[],[],Aeq,beq,zeros(1,H*N),ones(1,H*N),options);

%Gurobi solver
% params.outputflag = 0;          % Silence gurobi
% params.method     = 1;          % Use dual simplex method
% model_gurobi.A = sparse(Aeq);
% model_gurobi.obj = c;
% model_gurobi.sense = '=';
% model_gurobi.rhs = beq;
% result = gurobi(model_gurobi, params);
% x = result.x;
% Cmin = result.objval;

%sometimes some entries of matrix q can have very small negative values
%close to zero 
q = abs(reshape(x,N,H)');

end

