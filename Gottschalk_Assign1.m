%clear all past figures and data
clc;
clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #1 --Sample Space

%sample space
s = [2,3,4];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #2 --Outcomes
%number of trails
n = 30;
%generate random numbers for number of trials (n) either 2,3, or 4 in 
% row vector 'out'
out = randi([2,4],1,n);
%vector of 1 to number of trials conducted --- either 30,50,or 500
student = 1:n; 

%plot number of students/trials vs Outcomes
figure(1)
plot(student,out, 'r-o','Linewidth', 1.5)
xlabel('Number of Trials')
ylabel('Outcome')
title('Plot of Outcomes vs Trials')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #3 --Relative Frequency 
%intializing everything to zero -- f1,f2,f3 are vectors filled with zeros
f2=zeros(1,50);
f3=zeros(1,50);
f4=zeros(1,50);
n2=0;
n3=0;
n4=0;

%setting up vector and q number for the number of trials --changed from 50
%to 50 depending on question 3 or 4
q = 50;
trials = 1:q; 

%loop for amount of trials 
for i=1:q
    %set up random outcomes from 2,3,4 and put into vector outcome
    outcome = randi([2,4],1,i);

    %getting n value to calculate relative freq - find function gets matrix
    %of all the indexes of the number we are searching for in outcome vector
    %size function gets the size of the second number of the matrix -
    %therefore getting the sum of either 2,3,4 and puts sum in n2,n3, or n4
    n2= size(find(outcome==2),2);  
    n3= size(find(outcome==3),2);  
    n4= size(find(outcome==4),2); 

    %calc relative freq every loop by diving n2,n3,n4 by trial number (i)
    f2(i)=n2/i;
    f3(i)=n3/i;
    f4(i)=n4/i;
 end


%plot f2,f3,f4 in figure 2, added legend and x-y axias labels
figure(2)
plot(trials,f2,'r','LineWidth', 1)
hold on
plot(trials,f3,'g','LineWidth', 1)
hold on
plot(trials,f4,'b','LineWidth', 1)
hold on
xlabel('Number of Trials')
ylabel('Relative Frequency')
title('Plot of Relative Frequency vs Trials')
legend('2-C','3-B','4-A')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #4 --Relative Frequency
ff2=zeros(1,500);
ff3=zeros(1,500);
ff4=zeros(1,500);
nn2=0;
nn3=0;
nn4=0;

%setting up vector and q number for the number of trials --changed from 50
%to 500 depending on question 3 or 4
qq = 500;
trial = 1:qq; 

%loop for amount of trials 
for i=1:qq
    %set up random outcomes from 2,3,4 and put into vector outcome
    outcomes = randi([2,4],1,i);

    %getting n value to calculate relative freq - find function gets matrix
    %of all the indexes of the number we are searching for in outcome vector
    %size function gets the size of the second number of the matrix -
    %therefore getting the sum of either 2,3,4 and puts sum in nn2,nn3, or nn4
    nn2= size(find(outcomes==2),2);  
    nn3= size(find(outcomes==3),2);  
    nn4= size(find(outcomes==4),2); 

    %calc relative freq every loop by diving nn2,nn3,nn4 by trial number (i)
    ff2(i)=nn2/i;
    ff3(i)=nn3/i;
    ff4(i)=nn4/i;
 end


%plot ff2,ff3,ff4 in figure 2, added legend and x-y axias labels
figure(3)
plot(trial,ff2,'r','LineWidth', 1)
hold on
plot(trial,ff3,'g','LineWidth', 1)
hold on
plot(trial,ff4,'b','LineWidth', 1)
hold on
xlabel('Number of Trials')
ylabel('Relative Frequency')
title('Plot of Relative Frequency vs Trials')
legend('2-C','3-B','4-A')
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #5 --Pr[A]
%elementary vector of 1
elementary = 1:1; 

%dividing elementary vector divided by vector space to get probability of A
%given is that all 3 of sample space are equally likely
disp(length(elementary)/length(s));


