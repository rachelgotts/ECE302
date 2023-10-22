clc; close all; clear all;

num =[0 2];
dem = [1 25] ;

tr = 20*tf(num, dem)


% %%% feedback lead control
% num1 = [1 20];
% dem2 = [1 200];
% gs = 10*tf(num1,dem2)
% tr = tr*gs


figure(1)
subplot(2,1,1)
margin(tr)
grid on;
% subplot(2,1,2)
% nyquist(tr,'r')
% grid on;