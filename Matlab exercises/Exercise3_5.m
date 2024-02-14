%% Unconstrained optimization - Newton method with line search- Exercise 3.5

clear; close all;  

%% data 
alpha=0.1;
gamma=0.9;
tbar =1;
x0 = [0;0];
tolerance = 10^(-3) ;


x = x0 ;
X=[Inf,Inf,Inf,Inf,Inf];
for ITER=0:100
    [v, g, H] = f(x);
    X=[X;ITER,x',v,norm(g)];
    % stopping criterion
    if norm(g) < tolerance
        break
    end
    
    % search direction
    d = -inv(H)*g;
    t=tbar;
    while (f(x+t*d) > f(x)+alpha*t*d'*g)
        t=gamma*t;
    end
    % new point
    x = x + t*d; 

    ITER
    x
    v
    norm(g)
   
end

X


function [v, g, H] = f(x) 

v = 2*x(1)^4 + 3*x(2)^4 + 2*x(1)^2 + 4*x(2)^2 + x(1)*x(2) - 3*x(1) - 2*x(2)  ;

g = [ 8*x(1)^3 + 4*x(1) + x(2) - 3
      12*x(2)^3 + 8*x(2) + x(1) - 2];

H = [ 24*x(1)^2+4   1
          1        36*x(2)^2+8];

%% controllo calcolo gradiente ed hessiana
%syms x1 x2;
%v_symbolic = 2*x1^4 + 3*x2^4 + 2*x1^2 + 4*x2^2 + x1*x2 - 3*x1 - 2*x2;
    
% Calcolo del gradiente
%g_symbolic = gradient(v_symbolic, [x1, x2])
    
% Calcolo della matrice Hessiana
%H_symbolic = hessian(v_symbolic, [x1, x2])         

end
