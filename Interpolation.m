close all; clear all; clc;

% initial data
X = [2015; 2016; 2017];
Y = [31475; 30623; 28696];

% obtaining matrix A
A = [X(1)^0,X(1)^1,X(1)^2;
     X(2)^0,X(2)^1,X(2)^2;
     X(3)^0,X(3)^1,X(3)^2;]

% finding coeficients
c = A\Y;

% obtaining points 
XX = 2015:0.01:2017; %dense mesh for plot curve
YY = c(1)*XX.^0 + c(2)*XX.^1 + c(3)*XX.^2;

% ploting results;
figure(1); grid on; hold on;
plot(X, Y, '*b'); % ploting initial points
plot(XX, YY, '-r'); % ploting interpolating curve
 
 