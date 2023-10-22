%Gottschalk, Rachel ECE 302: Assignment #4
close all;
clear all;
clc;


%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%
% A - CDF
x = 0:((2*pi)/100):(2*pi); % x values for part 1 A
a = 0; % lower bound
b = 2*pi; % upper bound
y1= cdf1(x,a,b); % sends variables to function

%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%
y=0:0.03:0.98; % x values for part 2 A
x2=cdf2(y); % sends variables to function

%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%
x=0:0.02:10; % x values for part 3 A
x3=cdf3(x); % sends variables to function

%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%
function y = cdf1(x,a,b) 
    for i = 1:length(x) % finds length of vector x and loops through the length
        if x(i)<=0 % if x value is less than 0 then y=0
            y(i)=0;
        elseif x(i)>=a && x(i)<b % if x between 0 and 2pi then y is calulated
            y(i)=((x(i)-a)/(b-a));
        else % if x if greater than or equal to 2pi then y=1
            y(i)=1;
        end    
    end

    y1 = rand(100,1); % creates 100 random samples 
    x1=2*pi*y1; % puts random samples into inverse function to calculate x
    sample=1:100; % vector of 1 to 100
    
    %plots figure 1 with Uniform CDF and Inverse Function
    figure(1);
    subplot(2,1,1)
    plot(x,y,LineWidth=1.5,Color='cyan'); grid on;
    title('Uniform CDF - Part 1 (A)')
    xlabel("Theta")
    ylabel("F(Theta)")
    subplot(2,1,2)
    scatter(sample,x1,"cyan")
    title("100 Samples - Part 1 (B)")
    xlabel("Samples (n)")
    ylabel("Theta")
end

function x2 = cdf2(x)
    for i = 1:length(x) %finds length of vector x and loops through the length
        y(i) = ((x(i))^2+x(i))/2; % calculates the CDF of RV X
    end

    sample=1:100; % creates 100 random samples vector
    x2=sqrt(2*(rand(100,1))+1/4)-1/2; % puts random samples into inverse function to calculate x

    %plots figure 2 with Uniform CDF and Inverse Function
    figure(2);
    subplot(2,1,1)
    plot(x,y,LineWidth=1.5,Color='magenta'); grid on;
    title('CDF - Part 2 (A)')
    xlabel("X")
    ylabel("F(X)")
    subplot(2,1,2)
    scatter(sample,x2,"magenta")
    title("100 Samples - Part 2 (B)")
    xlabel("Samples (n)")
    ylabel("Random Variable (X)")

end
 
function x3 = cdf3(x)
    y=expcdf(x); %built in matlan function that calculates the exponential cdf of vector
    
    sample=1:100; % creates 100 random samples vector
    x3=-(log(1-rand(100,1))); % puts random samples into inverse function to calculate x

    %plots figure 1 with Expoential CDF and Inverse Function
    figure(3);
    subplot(2,1,1)
    plot(x,y,LineWidth=1.5,Color='green'); grid on;
    title('CDF - Part 3 (A)')
    xlabel("X")
    ylabel("F(X)")
    subplot(2,1,2)
    scatter(sample,x3,"green")
    title("100 Samples - Part 3 (B)")
    xlabel("Samples (n)")
    ylabel("Random Variable - X (in terms of 1/lamda)") %note lamda is just variable with no given values, so write RV X in terms of lamda

end

