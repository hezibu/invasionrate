function pst = count_p(t,params)
%This function calculates the value p from Solow and Costello (2004)
%It uses matrix coding for efficiency

S = repmat([1:t],t,1);
thing = 1-count_pi(S,S',params);
thing(t,:)=1; 
up = ones(size(thing)); upperones = triu(up,1); %upper triangular
thing2 = tril(thing) + upperones;
product = prod(thing2);

pst = product.*count_pi(1:t,t,params);