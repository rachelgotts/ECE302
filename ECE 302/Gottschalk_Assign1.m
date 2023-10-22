%clear all past figures and data
clc
clear all
close all 

%sample space
s = [2,3,4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #2 --Outcomes
%number of trails
n = 50;
%generate random numbers for number of trials (n) either 2,3, or 4 in 
% row vector 'out'
out = 1+randi(3,1,n);
%vector of 1 to number of trials conducted --- either 30,50,or 500
student = 1:n; 

%plot number of students/trials vs Outcomes
figure(1)
plot(student,out, 'r-o','Linewidth', 1.5)
xlabel('Number of Trials')
ylabel('Outcome')
title('Plot of Outcomes vs Trials')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #3&4 --Relative Frequency 
%intializing everything to zero -- f1,f2,f3 are vectors filled with zeros
f2=zeros(1,50);
f3=zeros(1,50);
f4=zeros(1,50);
n2=0;
n3=0;
n4=0;

%loop for 50 trials 
for i=1:50
    %if statements check the numbers in 'out' array and add up the amount 
    % of each to get N1(i),N2(i),N3(i) values
    if out(i)==2
        n2=n2+1;
    elseif out(i)==3
        n3=n3+1;
    else
        n4=n4+1;
    %calc relative freq every loop by diving n2,n3,n4 by trial number
    f2(i)=n2/i;
    f3(i)=n3/i;
    f4(i)=n4/i;
    end

end

%plot f2,f3,f4 in figure 2, added legend and x-y axias labels
figure(2)
plot(student,f2,'r')
hold on
plot(student,f3,'g')
hold on
plot(student,f4,'b')
hold on
xlabel('Number of Trials')
ylabel('Relative Frequency')
title('Plot of Relative Frequency vs Trials')
legend('2-C','3-B','4-A')
grid on

