function Gottschalk_Assign7
% Sample MATLAB code template for implementing KF for 1D vehicle navigation.
% This sample code may serve as a guide but feel free to implement using
% your own coding style. 
clear;
close all;
clc;
% ---------------------
% DEFINE SIMULATION PARAMETERS
% ---------------------
% For example: m, tf, T, tSteps, N=1, sigmaM, R , q_0, q_0HatPlus, u,
% P_0Plus, Q, R= sigmaM^2, F, G, H 
m = 1; % mass
tf = 60; % final time 
T = 0.5; % fixed time interval
time = 0:T:tf; 
sigmaM = 0.5; % sigma
tSteps = size(time,2);
u = 0.1; % force
R = sigmaM^2; % sigma ^2
N=1;

x01 = 2; % true postion at t=0
x0e = 3; % estimated position at t=0
v01 = 4; % true velcocity at ms^-1
v0e = 6; % estimated velcoity

q_0 = [2 4]';
q_0HatPlus = [3 6]';
P0 = diag([(x01-x0e)^2,(v01-v0e)^2]); 

F = [1 T; 0 1];
G = [0 T/m]';
H = [0 1]; 
Q = 10^(-2)*[0.2 0.01; 0.01 0.1];

% --------------
% RUN SIMULATION
% --------------
% Update the state equations

for k = 1:tSteps   % Main timing loop  (you will need to write code for this loop)
    
    % Run the true (without noise) vehicle model
    q(:,k) = F*q_0 + G*u;

    % %%%%%%%%%%% Run the KF estimator
        
    % Generate process noise, zeta 
    [U,D,V]=svd(Q);
    zeta = V*(D^(0.5))*(randn(2,1)); % error

    % Prediction step (predicted state and covariance)
    q_hat_predict = F*q_0HatPlus + G*u + zeta; 
    P_minus_predict = F*P0*F' + Q; 
    
    % Measurement with noise 
    zi = randn(1,1);  
    z = H*q(:,k)+zi; 
    y = H*q(:,k); 

    % Compute the Kalman filter gain, K_kPlus1
    K = P_minus_predict*H'*inv(H*P_minus_predict*H'+R); %kalman gain
    
    % Update step (updated/corrected state and covariance)
    q_hat(:,k) = q_hat_predict+K*(z-H*q_hat_predict);
    P(:,:,k) = (eye(2) - K*H)*P_minus_predict;

    % Make sure to store updated state vector and covariance matrices.
    q_0 = q(:,k);
    q_0HatPlus = q_hat(:,k);
    P0 = P(:,:,k);
                                
end
% -----------
% PLOT OUTPUT
% -----------
tOut = T*(1:tSteps);
% Plot the estimation error with covariance bounds
figure;
clf;
subplot(2,1,1);
hold on;
plot(tOut,q(1,:)-q_hat(1,:),'b');
sigma(1,:) = P(1,1,:);
plot(tOut,3*sqrt(sigma(1,:)),'g');
plot(tOut,-3*sqrt(sigma(1,:)),'g');
hold off;
box on;
ylabel('e_x [m] \pm 3\sigma_1');
subplot(2,1,2);
hold on;
plot(tOut,q(2,:)-q_hat(2,:),'b');
sigma(2,:) = P(2,2,:);
plot(tOut,3*sqrt(sigma(2,:)),'g');
plot(tOut,-3*sqrt(sigma(2,:)),'g');
hold off;
box on;
ylabel('e_v [m/s] \pm 3\sigma_2');
xlabel('Time [s]');
% savefilename = 'error';
% exportgraphics(gca,[savefilename,'.pdf'],'Resolution',300);
%saveas(gcf, savefilename, 'fig');
%print('-depsc2', '-r300', [savefilename, '.eps']);
% Create a movie of the simulation
figure;
clf reset;
xmax = max(q(1,:));
xmin = min(q(1,:));
ymax = max(q(2,:));
ymin = min(q(2,:));
j = 0;
for k = 1:5:tSteps
    clf;
    box on;
    axis([xmin-5 xmax+5 ymin-5 ymax+5]);
    axis equal;
    axis manual;
    plot(q(1,k),q(2,k),'>');
    hold on;
    
    desired = plot(q(1,1:k),q(2,1:k),'b-');
    estimated = plot(q_hat(1,1:k),q_hat(2,1:k),'g.-');
    
    
    [U,S,V] = svd(P(:,:,k));
    ellipse(10*S(1,1),10*S(2,2),atan2(U(1,2),U(1,1)),q_hat(1,k),q_hat(2,k),'r');
    hold off;
    xlabel('x - position [m]');
    ylabel('v - velocity [m/s]');
    j = j + 1;
    %drawnow;
    
    M1(j) = getframe;
end
legend([desired estimated], 'Ref. trajectory', 'Estimated trajectory','Location', 'Best');
grid on
% savefilename = 'OUT/trajectorySim';
% exportgraphics(gcf,[savefilename,'.pdf'],'Resolution',300);
% saveas(gcf, savefilename, 'fig');
% print('-depsc2', '-r300', [savefilename, '.eps']);
disp('... simulation complete. ');
function h=ellipse(ra,rb,ang,x0,y0,C,Nb)
% Ellipse adds ellipses to the current plot
%
% ELLIPSE(ra,rb,ang,x0,y0) adds an ellipse with semimajor axis of ra,
% a semimajor axis of radius rb, a semimajor axis of ang, centered at
% the point x0,y0.
%
% The length of ra, rb, and ang should be the same. 
% If ra is a vector of length L and x0,y0 scalars, L ellipses
% are added at point x0,y0.
% If ra is a scalar and x0,y0 vectors of length M, M ellipse are with the same 
% radii are added at the points x0,y0.
% If ra, x0, y0 are vectors of the same length L=M, M ellipses are added.
% If ra is a vector of length L and x0, y0 are  vectors of length
% M~=L, L*M ellipses are added, at each point x0,y0, L ellipses of radius ra.
%
% ELLIPSE(ra,rb,ang,x0,y0,C)
% adds ellipses of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% ELLIPSE(ra,rb,ang,x0,y0,C,Nb), Nb specifies the number of points
% used to draw the ellipse. The default value is 300. Nb may be used
% for each ellipse individually.
%
% h=ELLIPSE(...) returns the handles to the ellipses.
%
% as a sample of how ellipse works, the following produces a red ellipse
% tipped up at a 45 deg axis from the x axis
% ellipse(1,2,pi/8,1,1,'r')
%
% note that if ra=rb, ELLIPSE plots a circle
%
% written by D.G. Long, Brigham Young University, based on the
% CIRCLES.m original 
% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch
% Check the number of input arguments 
if nargin<1,
  ra=[];
end
if nargin<2,
  rb=[];
end
if nargin<3,
  ang=[];
end
%if nargin==1,
%  error('Not enough arguments');
%end
if nargin<5,
  x0=[];
  y0=[];
end
 
if nargin<6,
  C=[];
end
if nargin<7,
  Nb=[];
end
% set up the default values
if isempty(ra),ra=1;end
if isempty(rb),rb=1;end
if isempty(ang),ang=0;end
if isempty(x0),x0=0;end
if isempty(y0),y0=0;end
if isempty(Nb),Nb=300;end
if isempty(C),C=get(gca,'colororder');end
% work on the variable sizes
x0=x0(:);
y0=y0(:);
ra=ra(:);
rb=rb(:);
ang=ang(:);
Nb=Nb(:);
if isstr(C),C=C(:);end
if length(ra)~=length(rb),
  error('length(ra)~=length(rb)');
end
if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end
% how many inscribed elllipses are plotted
if length(ra)~=length(x0)
  maxk=length(ra)*length(x0);
else
  maxk=length(ra);
end
% drawing loop
for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    radm=ra(k);
    radn=rb(k);
    if length(ang)==1
      an=ang;
    else
      an=ang(k);
    end
  elseif length(ra)==1
    xpos=x0(k);
    ypos=y0(k);
    radm=ra;
    radn=rb;
    an=ang;
  elseif length(x0)==length(ra)
    xpos=x0(k);
    ypos=y0(k);
    radm=ra(k);
    radn=rb(k);
    an=ang(k)
  else
    rada=ra(fix((k-1)/size(x0,1))+1);
    radb=rb(fix((k-1)/size(x0,1))+1);
    an=ang(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end
  co=cos(an);
  si=sin(an);
  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
%  x=radm*cos(the)*co-si*radn*sin(the)+xpos;
%  y=radm*cos(the)*si+co*radn*sin(the)+ypos;
  h(k)=line(radm*cos(the)*co-si*radn*sin(the)+xpos,radm*cos(the)*si+co*radn*sin(the)+ypos);
  set(h(k),'color',C(rem(k-1,size(C,1))+1,:));
end