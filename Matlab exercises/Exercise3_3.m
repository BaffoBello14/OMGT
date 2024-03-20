%% Unconstrained optimization -- Exercise 3.3

clear; close all; clc; 


%% The problem: min f(x(1),x(2))= x(1)^4 + x(2)^4 - 2*x(1)^2 + 4*x(1)*x(2)-2*x(2)^2

%% Data

alpha = 0.1;
gamma = 0.9;
tbar = 1;
x0 = [ 0 ; 0; 0];
tolerance = 10^(-3) ;

%% Method: gradient method with inexact line search


X=[Inf,Inf,Inf,Inf,Inf];

ITER = 0 ;
x = x0 ;

while true
    [v, g] = f(x);
    
    X=[X;ITER,x(1),x(2),v,norm(g)];
    
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    % search direction
    d = -g;
    
    % Armijo inexact line search
    t = tbar ;
    while f(x+t*d) > v + alpha*g'*d*t
        t = gamma*t ;
    end
        
    % new point
    x = x + t*d ;
    ITER = ITER + 1 ;


    
end
disp('optimal solution')
x
v
norm(g)
ITER
function [v, g] = f(x) 

v = exp(-x(1)-x(2)-x(3))+ x(1)^2 + 3*x(2)^2 + x(3)^2 + x(1)*x(2)  - x(2)*x(3) + x(1) - 3*x(3);

%% controllo calcolo gradiente ed hessiana
syms x1 x2 x3;
v_symbolic = exp(-x1-x2-x3)+ x1^2 + 3*x2^2 + x3^2 + x1*x2 - x2*x3 + x1 - 3*x3;

    
% Calcolo del gradiente
g_symbolic = gradient(v_symbolic, [x1, x2, x3])
    
% Calcolo della matrice Hessiana
%H_symbolic = hessian(v_symbolic, [x1, x2]) 

g = [ 2*x(1) + x(2) - exp(- x(1) - x(2) - x(3)) + 1
x(1) + 6*x(2) - x(3) - exp(- x(1) - x(2) - x(3))
 2*x(3) - x(2) - exp(- x(1) - x(2) - x(3)) - 3];

end



