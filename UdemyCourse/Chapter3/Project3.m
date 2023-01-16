% Function to be implemented in simulink:
% A dx/t + Ti int(x) dt = F

% bring all the signals to the left, let it remain just the derivative
% dx/t  = (F - + Ti int(x) dt) / A

%  Here, we define the parameters

F = 2;
A = 3;
Ti = 5;

sim('Project3Simulink.slx')