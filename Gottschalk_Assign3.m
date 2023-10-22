%Gottschalk, Rachel ECE 302: Assignment #3
close all;
clear all;
clc;


n = 100; % max of sample number
xs = 1:n; % sample number array

%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig1=1; % figure 1 variable
m = 1; % mean 
o = 0.3; % standard deviation
x1 = (m + o * randn(n,1))'; % generates random sample measurements x
x11=sort(x1); % sorts x1 from smaller to larger values

a(m,o,x11,fig1); % passes mean, standard deviation, sorted x, and fig1 to a function
b(x1,xs,o,m,fig1); % passes unsorted x, sample number array, standard deviation, mean and fig1 to b function


%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig2=2; % figure 2 variable
m = 3; % mean 
o = 0.8; % standard deviation
x2 = (m + o * randn(n,1))'; % generates random sample measurements x
x22 = sort(x2); % sorts x1 from smaller to larger values

a(m,o,x22,fig2); % passes mean, standard deviation, sorted x, and fig1 to a function
b(x2,xs,o,m,fig2); % passes unsorted x, sample number array, standard deviation, mean and fig1 to b function




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fx = a(m,o,x,fig) 
    for i = 1:length(x) % finds length of x and makes array from 1 to that length and iterates through that array
        fx(i) = (1/(o*sqrt(2*pi)))*exp(-0.5*(((x(i)-m)^2)/(o^2))); % gaussian fx(x) eq 
    end
    %plots figure 1 with the first subplot 
    figure(fig) 
    subplot(2,1,1)
    plot(x,fx, 'r',"LineWidth", 1.5)
    title('PDF of X')
    xlabel('X - Sample Measurements')
    ylabel('PDF')
    grid on
    hold on
end


function y = b(x,xs,o,m,fig)
    % makes confidence interval lines 
    y=m+o; 
    y1=m-o;
    y2=m+2*o;
    y3=m-2*o;
    y4=m+3*o;
    y5=m-(3*o);

    % plots scatter of N vs X and has confidence intervals on same subplot
    figure(fig)
    subplot(2,1,2)
    scatter(xs,x, 'b')
    hold on
    yline(y)
    hold on
    yline(y1)
    hold on 
    yline(y2,LineWidth=1)
    hold on 
    yline(y3',LineWidth=1)
    hold on 
    yline(y4,LineWidth=1.5)
    hold on 
    yline(y5,LineWidth=1.5)
    hold on
    title('Sample Number (N) vs Sample Measurements (X)')
    xlabel('N - Sample Numbers')
    ylabel('X - Sample Measurements')
    legend('N vs X','m+o','m-o','m+2o','m-2o','m+3o','m-3o')
    grid on
 
end




