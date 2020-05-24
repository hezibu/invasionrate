function pi = count_pi(s,t,params)
%This function calculates the variable pi from Solow and Costello (2004)
%disp('new pi')
pi0 = params(3);
pi1 = params(4);
pi2 = params(5);

num = exp(pi0 + pi1*t + pi2*exp(t-s));
f = find(isinf(num)); %if the number is infinity
pi = num./(1+num);
pi(f) = 1;

