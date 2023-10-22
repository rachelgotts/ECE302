%Gottschalk, Rachel ECE 302: Assignment #2
close all;
clear all;
clc;

%%%%%%%%%%%%%%%% 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RV X = Number of pounds rounded
% RV Y = Charge in Cents for sending 1 package --- use function y=g(x)
% $1 = 1pound, $0.90 = 2pound, $0.80 = 3pounds, $0.70 = 4pounds, 
% $0.60 = 5pounds, $5.00 6pounds<=X<=10pounds, X>10 = will not accept4

x = 10*rand; % select rand numbers for x 
disp(g(1)); % pass 5 through the function to show that it is grabbing the x values and then plugging them into g(x)

%%%%%%%%%%%%%%%% 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sx = [1 2 3 4]; % initialize range of x
Sy = [100 190 270 340]; % intialize range of y

% loop going each of the Sx and Sy to calculate the PMF of them
for i = 1:length(Sx)
    %%%Px(x) x-g(x)
    % Discrete Uniform RV - X = 1,2,3,4 - so use 1/(l-k+1)
    px(i) = 1/(length(Sx)-Sx(1)+1);
    % Discrete Uniform RV - Y = 1, 1.9, 2.7, 3.4 - so use 1/(l-k+1)
    py(i) = 1/(length(Sy)-1+1);
end

% plotting the PMG of Px and Py
figure(1)
subplot(2,1,1)
stem(Sx,px,'b')
title('PMF of X')
xlabel('RV (X)')
ylabel('PMF X')
hold on;
grid on;
subplot(2,1,2)
stem(Sy,py,'r')
title('PMF of Y')
xlabel('RV (Y)')
ylabel('PMF Y')
grid on;
hold on;


Ey=0;
% loop through Sy to calulate expected values
for i = 1:length(Sy)
    Ey = Ey + Sy(i)*py(i);
end

%%%%%%%%%%%%%%%% 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sx3 = 1:8; % intializing range of RV x to be from 1 to 8
Sy3 = [100 190 270 340 400 500]; % initializing range of RV y 

py=(zeros(size(Sy3))); % setting values (number is the amount of Sy3) of py to all zeros and
Ey=0; % initializing expected value of y to 0

py5=0;
% looping through all of the values in the range of Sx3
for i = 1:length(Sx3)
    if 0<i && i<5 % checking to see if value of Sx3 is greater than 0 and less than 5
        px(i)= 0.15; 
        py(i) = px(i);
    elseif 5<i && i<9 % checks if i is greater than 5 and less than 9
        px(i)=0.1;
        py5=py5+0.1;
    else % if i = 5
        px(i) = 0.1;
        py(i) = px(i);
    end
   py(6)=py5;
end

Ey = (Sy3(1) + Sy3(2) + Sy3(3) + Sy3(4))*px(1) + Sy3(5)*px(6)+Sy3(6)*py(6); % expected value of y

figure(2)
subplot(2,2,1)
stem(Sx3,px,'b')
title('PMF of X')
xlabel('RV (X)')
ylabel('PMF')
hold on;
grid on;
subplot(2,2,2)
stem(Sy3,py,'r')
title('PMF of Y')
xlabel('RV (Y)')
ylabel('PMF')
hold on;
grid on;

%%%%%%%%%%%%%%%% 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(shipweight8(10)); 

%%%%%%%%%%%%%%%% 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(avg(10));
disp(avg(100));
disp(avg(1000));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% 5 - Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=avg(m)
    s = cumsum([100 90 80 70 60 50]); % array of cumulative sum starting at the beginning of the first array 
    y2 = [s 500 500 500]; % places s in array y2 
    y1=sum(y2); % sums the y2 vector 
    y=y1/m; % divides y1 by m to get the average
end

%%%%%%%%%%%%%%%% 4 - Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,x] = shipweight8(m)
    rx = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 6 6 7 7 8 8]; % takes in account that 3/20 is the same probility and then we normalize
    x=zeros(1,m); % puts zeros in x for all m
    y=zeros(1,m); % puts zeros in y for all m

    for n = 1:m
        x(n)=rx(randi(length(rx))); % picks randi number of the length of rx and grabs that correspodning number in rx
        y(n)=g(x(n)); % sends value of x(n) to g(x) function and returns corresponding y value
    end

% in class:
%     sx = 1:8;
%     p = [0.15 0.15 0.15 .15 0.1 0.1 0.1 0.1];
%     
%     c =cumsum([0,p(:).']);
%     c = c/c(end); % take last value of index and divide entire cummlative sum
% 
%     randNo = rand(1,m);
%     [N,i]=histc(randNo,c);
%     x=sx(i); % map to v values
end


%%%%%%%%%%%%%%%% 1 - Function %%%%%%%%%%%%%%%%%%

function y = g(x)
    x = ceil(x);

    if x > 10
        y = NaN; 
    end

    switch x 
        case 1
            y = 100;
        case 2
            y = 190;
        case 3
            y = 270;
        case 4
            y = 340;
        case 5
            y = 400;
        case 5<x && x<11
            y= 500;
        otherwise
            y= 0;
            disp("Pounds not Accepted"); 
      
    end
end




