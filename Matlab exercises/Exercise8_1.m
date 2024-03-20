%% Constrained optimization -- Exercise 8.1

clear; close all; clc;

%% data 

global Q c A b eps;

%% data

Q = [ 1 0 ; 0 2 ] ;
c = [ -3 ; -4 ] ;
A = [2 1; -1 -1; -1 0 ];
b = [ 4;-1;0];

tau = 0.1 ;
eps0 = 5 ;
tolerance = 1e-6 ;

%% penalty method


eps = eps0;
x = [4;0];
iter = 0;
SOL=[];

while true
    [x,pval] = fminunc(@p_eps,x);      
    infeas = max(A*x-b);
    
    SOL=[SOL;iter,eps,x',infeas,pval];
    if infeas < tolerance
        break
    else
        eps = tau*eps;
        iter = iter + 1 ;
    end
end



fprintf('\t iter \t eps \t x(1) \t x(2) \t max(Ax-b) \t pval \n');

SOL
%% penalized function

function v= p_eps(x) 

    global Q c A b eps;

    %v = 0.5*x'*Q*x + c'*x ;

    % se la funzione non è scrivibile in forma matriciale come sopra (in
    % forma quadratica), allora scriverla semplicemente nella forma
    % sottostante,cioè v = funzione

    v = x(1)^2 - log(x(1)+x(2));
  

    for i = 1 : size(A,1)
        v = v + (1/eps)*(max(0,A(i,:)*x-b(i)))^2 ;
    end

end