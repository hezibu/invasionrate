%This code performs the main calculations for Solow and Costello (2004) 

clear all
close all

global num_discov; %will need to be called by some function files
load NumDis.txt %load the text file of the number of discoveries by year
T=[1851:1995]; %the time period over which discoveries were made
num_discov = NumDis';
disp('SF data set')
options = optimset('TolFun',.01,'TolX',.01); %set the tolerances 

guess =    [-1.1106;    0.0135;   -1.4534; 0; 0]; %an initial guess 
constr = 99*ones(size(guess)); %to set a constraint make some of these different than 99

[vec1 val1] = fminsearch('count_log_like',guess,options,constr) %Matlab routine for conducting MLE
C1 = count_lambda(vec1,length(num_discov)); %Calculates the mean of Y

%Create the plot 
plot(T,cumsum(num_discov),'k-',T,cumsum(C1),'k--')
legend('Discoveries','Unrestricted')
xlabel('Year')
ylabel('Cumulative Discovery')

