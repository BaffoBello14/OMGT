%% Classification problems - Exercise 5.5

close all; clear; clc;

%% data

A=[ 0.0113    0.2713
    0.9018   -0.1121
    0.2624   -0.2899
    0.3049    0.2100
    -0.2255   -0.7156
    -0.9497   -0.1578
    -0.6318    0.4516
    -0.2593    0.6831
    0.4685    0.1421
    -0.4694    0.8492
    -0.5525   -0.2529
    -0.8250    0.2802
    0.4463   -0.3051
    0.3212   -0.2323
    0.2547   -0.9567
    0.4917    0.6262
    -0.2334    0.2346
    0.1510    0.0601
    -0.4499   -0.5027
    -0.0967   -0.5446];

B=[ 1.2178    1.9444
    -1.8800    0.1427
    -1.6517    1.2084
    1.9566   -1.7322
    1.7576   -1.9273
    0.7354    1.1349
    0.1366    1.5414
    1.5960    0.5038
    -1.4485   -1.1288
    -1.2714   -1.8327
    -1.5722    0.4658
    1.7586   -0.5822
    -0.3575    1.9374
    1.7823    0.7066
    1.9532    1.0673
    -1.0233   -0.8180
    0.8021    0.3341
    0.0473   -1.6696
    0.8783    1.9846
    -0.5819    1.8850];

nA = size(A,1);
nB = size(B,1);

% training points
T = [A ; B]; 

y = [ones(nA,1) ; -ones(nB,1)]; % labels
l = length(y);

%% Nonlinear SVM

% parameter
C = 1 ;

% Gaussian kernel
gamma = 1 ;
K = zeros(l,l);
for i = 1 : l
    for j = 1 : l
        K(i,j) = exp(-gamma*norm(T(i,:)-T(j,:))^2);
    end
end

% Linear kernel
% K = T * T';

% Polinomial kernel
% p = 2; % Imposta il grado del polinomio
% K = (T * T' + 1).^p;

% Tanh kernel
% beta = 1;
% gamma = 0;
% K = tanh(beta * (T * T') + gamma);

% define the problem
Q = zeros(l,l);
for i = 1 : l
    for j = 1 : l
        Q(i,j) = y(i)*y(j)*K(i,j) ;
    end
end

% solve the problem
[la,ov] = quadprog(Q,-ones(l,1),[],[],y',0,zeros(l,1),C*ones(l,1));

% compute b
ind = find((la > 1e-3) & (la < C-1e-3));
i = ind(1);
b = 1/y(i) ;
for j = 1 : l
    b = b - la(j)*y(j)*K(i,j);
end

% compute w
w = zeros(size(T(1,:)));
for i = 1:l
    w = w + la(i) * y(i) * T(i,:);
end


%% Write explicitly the vector of the optimal solution of the dual problem
% Il vettore dei moltiplicatori di Lagrange ottimali è 'la'
% peso ottimale 'w'
% offset ottimale 'b'

la
w
b

%% Find the misclassified points of the data sets A and B by means of the dual solution
% Calcola la funzione decisionale per tutti i punti
f = zeros(l, 1);
for i = 1:l
    f(i) = sum(la .* y .* K(:, i)) + b;
end

% Calcola gli errori di classificazione xi per tutti i punti
xi = max(0, 1 - y .* f)

% Trova i punti misclassificati
misclassified_indices = find(xi > 1)
misclassified_points = T(misclassified_indices, :)

% Compute the support vectors

supp = find(la > 10^(-3));
suppA = supp(supp <= nA);
suppB = supp(supp > nA);


%% (d) Classify the new point (0,1)
%% caso gaussiano
new_point = [0, 1];
s = 0;
for i = 1:l
    s = s + la(i) * y(i) * exp(-gamma * norm(T(i,:) - new_point)^2);
end
s = s + b

%% caso lineare
%s = 0;
%for i = 1:l
%    s = s + la(i) * y(i) * (new_point * T(i,:)');
%end
%s = s + b;

%% caso polinomiale
%s = 0;
%for i = 1:l
%    s = s + la(i) * y(i) * ((new_point * T(i,:)' + 1)^p);
%end
%s = s + b;

% caso tanh
%s = 0;
%for i = 1:l
%    s = s + la(i) * y(i) * tanh(beta * (new_point * T(i,:)') + gamma);
%end
%s = s + b;


if s > 0
    fprintf('Il nuovo punto appartiene alla classe A.\n');
else
    fprintf('Il nuovo punto appartiene alla classe B.\n');
end

%% plot the surface f(x)=0
for xx = -2 : 0.01 : 2
    for yy = -2 : 0.01 : 2
        s = 0;
        for i = 1 : l
            s = s + la(i)*y(i)*exp(-gamma*norm(T(i,:)-[xx yy])^2);
            %s = s + la(i)*y(i)*([xx yy] * T(i,:)'); lineare
            %s = s + la(i)*y(i)*((([xx yy] * T(i,:)') + 1)^p); polinomiale
            %s = s + la(i)*y(i)*tanh(beta * ([xx yy] * T(i,:)') + gamma); tanh
        end
        s = s + b;
        if (abs(s)< 10^(-2))
            plot(xx,yy,'g.');
            hold on
        end
    end
end
plot(A(:,1),A(:,2),'bo',B(:,1),B(:,2),'ro','Linewidth',5)
