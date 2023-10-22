%Gottschalk, Rachel ECE 302: Assignment #5
close all;
clear all;
clc;


%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%
m = 100; 
x = linspace(-3,3,m);
y = linspace(-3,3,m);

fxy=biPdf(x,y,0,0,1,1,0.6); % sends to function

%plots part 1 with mechc
figure(1) 
meshc(x,y,fxy)
title("Part 1: Standard Bivariate Gaussian PDF")
xlabel('x');
ylabel('y');
zlabel('Joint PDF f_x_,_y(x,y)');
grid on;


%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%
v = linspace(-3,3,m);
w = linspace(-3,3,m);
M=[v;w];
A =  1/sqrt(2)*[1 1; -1 1];   %    A =  [1/sqrt(2)  1/sqrt(2)]
                              %         [-1/sqrt(2) 1/sqrt(2)]
% calculates transform
t = A^-1*M;
x = t(1,:);               
y = t(2,:);

% sends to function and then divides by determinate 
fvw = biPdf(x,y,0,0,1,1,0.6);
fvw = fvw/det(A);

%plots figure 2
figure(2)
meshc(v,w,fvw)
title('Part 2: PDF of linear transformations')
xlabel('v');
ylabel('w');
zlabel('Joint PDF f_v_,_w(v,w)')


%%%%%%%%%%%%%%%%%%%% Part 3 %%%%%%%%%%%%%%%%%%%%%%
mx = 3;
ox = 5^2;
my = 3;
oy = 1^2;
p = 0.7;

% generates random numbers and then makes A matrix 
z = randn(2,m);
A = [sqrt(ox) 0; p*sqrt(oy) sqrt((1-p^2)*oy)];
M = [mx;my];
X = M + A*z;
x = X(1,:);
y = X(2,:);

% sorts x and y 
x1 = sort(x);
y1 = sort(y);

%sends to function to calculate buvariate PDF
fxy = bivariatePdf(x1,y1,mx,ox,my,oy,p);

% plots figure 3 with all parts on one figure
figure(3)
subplot(2,1,1)
scatter(x,y,'filled')
title('Part 3a: Generating Bivariate Samples Using Standard Gaussian PDF');
xlabel('x');
ylabel('y');
subplot(2,1,2)
meshc(x1,y1,fxy)
title('Part 3b: Generating Bivariate Samples Using Standard Gaussian PDF');
xlabel('x');
ylabel('y');
zlabel('Joint PDF f_x_,_y(x,y)')


%%%%%%%%%%%%%%%%%%%% Part 4 %%%%%%%%%%%%%%%%%%%%%%
% intial values
mx = 0;
sigx = 2;
ox=4;
mn = 0;
sign = 2;
on=4;
my = 0;
sigy=sqrt(ox+on);
oy=sigy^2;
p = 0;

% generates normal distrubuted numbers
x = normrnd(0,sigx,1,m);
y = normrnd(0,sigy,1,m);

% stores origional numbers
x0=x;
y0=y;

% sort x and y and make A matrix
x=sort(x);
y=sort(y);
A = [1 0; 1 1];

% makes function to calculate x and n matrixes
t=A\[x;y];
x = t(1,:);
n = t(2,:);

% sends to function
fxy = bivariatePdf(x,n,mx,ox,mn,on,p);
fxy = fxy/det(A);

% plots all parts of 4 on one figure
figure(4)
subplot(2,2,1)
meshc(x,y,fxy)
title('Part 4a: Estimating Signal from Noisy Measurements')
xlabel('x');
ylabel('y');
zlabel('Joint PDF f_x_,_y(x,y)')

% generated x and n
tran = A\[x0;y0];
x = tran(1,:);
n = tran(2, :);

% stores origin variables
x1 = x;
n1 = n;

y1 = y0;

% created estimated number
x_estimated = (ox/(ox+on))*y1;

% plots part 4b
vals = 1:1:m;
subplot(2,2,2)
plot(vals, x1,'g')
hold on;
plot(vals, y1,'b-')
hold on;
plot(vals,x_estimated,"k-")
hold on;
title('Part 4b: Samples of the Original, Measured, and Estimated Signals')
xlabel('Sample Number'); ylabel('Signal');
legend('Original Signal','Measured Signal','Estimated Signal');

% generated errors
error_M = abs(y1-x1);
error_E = abs(x_estimated - x1);

plots part 4c
subplot(2,1,2)
plot(vals,error_M,"r")
hold on;
plot(vals, error_E, "b")
title("Part 4c: Errors")
xlabel('Sample Number'); ylabel('Error')
legend("Measured Error", "Estimated Error")


% part d - calulates sum and prints them in command window
s1 = sum((y1-x1).^2)
s2 = sum((x_estimated-x1).^2)



%%%%%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function fxy = biPdf(x,y,mx,my,ox,oy,p) 

    % calculates the bivar PDF with 1st formula
    for i = 1:length(x)
        for j = 1:length(y)
            fxy(i,j) = 1/((ox*ox)*sqrt(1-p^2)*sqrt((2*pi)^2))*exp(-.5*(((x(i)-mx)/ox)^2-(2*p*(x(i)-mx)/ox)*((y(j)-my)/oy)+(y(j)-my/oy)^2)/(1-p^2));
        end
    end

end 


function f_xy = bivariatePdf(x,y,mu1,var1,mu2,var2,corr_coeff)
    cov = [var1,  corr_coeff*sqrt(var1)*sqrt(var2); % Calculate covariance matrix
           corr_coeff*sqrt(var1)*sqrt(var2)   var2];
  
    mu = [mu1;mu2]; % Mean matrix

    f_xy = zeros(length(x),length(y));

    for i=1:length(x) % Nested for loop to calculate bivariate guassian pdf
        for j=1:length(y)
            xVec = [x(i); y(j)];
            f_xy(i,j) = (1/(2*pi*sqrt(det(cov))))*exp(-0.5*(xVec-mu)'*cov^-1*(xVec-mu));
        end
    end
end

















