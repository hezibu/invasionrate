function [lambda] = count_lambda(params,N)
%This function calculates lambda from
%Solow and Costello, 2004.
%params is a vector of parameters

for t=1:N
    S=1:t;
    Am = count_m(S,params); 
    Ap = count_p(t,params);
    lambda(t) = Am*Ap';
end
 
