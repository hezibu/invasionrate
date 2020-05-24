function [LL, lambda] = count_log_like(params,restrict)
%This function file calculates the log likelihood function for Solow and
%Costello (2004).  It takes into account any possible restrictions (See
%below)

%params is a vector of parameters
%restrict is a vector (same size as params) that places restrictions on the parameters.  If restrict(i)=99, then there
%is no restriction for the ith parameter.  If restrict(i)=0 (for example)
%then the restriction is exactly that.

global num_discov;

f=find(restrict~=99);
g=find(restrict==99);
params(g) = params;
params(f) = restrict(f);

%
for t=1:length(num_discov)
    S=1:t;
    Am = count_m(S,params); 
    Ap = count_p(t,params);
    lambda(t) = Am*Ap';
    summand2(t) = num_discov(t)*log(lambda(t)) - lambda(t);
end

LL = -sum(summand2)
