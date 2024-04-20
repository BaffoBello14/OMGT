function bimatrix_games()

    clc, clear all

    P1 = [
    4 3 2
    2 1 5
    ];

    P2 = [
    5 4 3
    7 2 6
    ];

    stable = false;
    
    salva1 = P1;
    salva2 = P2;

    while ~stable
        [P1_new, P2_new] = strictly_dominated_strategies(P1, P2);
        
        stable = isequal(P1, P1_new) && isequal(P2, P2_new);

        P1 = P1_new;
        P2 = P2_new;
    end

    P1 = salva1;
    P2 = salva2;
    
    disp('Reduced matrix P1:');
    disp(P1_new);
    
    disp('Reduced matrix P2:');
    disp(P2_new);
    
    pureNashEquilibria(P1, P2);

    disp("Mixed for P1:")
    mixedNashEquilibria(P1_new, 0);

    disp("Mixed for P2:")
    mixedNashEquilibria(P2_new, 1);

    disp("Mixed with KKT:")
    for i = 1:10 % Viene fatto più volte così da trovare tutti i casi
        mixedKKT(P1, P2);
    end
end

function [P1_new, P2_new] = strictly_dominated_strategies(P1, P2)
    % Inizializza le matrici P1_new e P2_new con le stesse dimensioni di P1 e P2
    P1_new = P1;
    P2_new = P2;

    % Elimina le righe in base alla condizione specificata
    for i = 1:size(P1, 1)
        if any(all(P1 < P1(i, :), 2))
            disp(['Strategy ', num2str(i), ' for Player 1 is strictly dominated by strategy ', num2str(find(all(P1 < P1(i, :), 2))), '.']);
            P1_new(i, :) = [];
            P2_new(i, :) = [];
            % Dopo l'eliminazione di una riga, ricomincia il loop
            % dalla prima riga
            break;
        end
    end

    % Elimina le colonne in base alla condizione specificata
    for i = 1:size(P2, 2)
        if any(all(P2 < P2(:, i), 1))
            disp(['Strategy ', num2str(i), ' for Player 2 is strictly dominated by strategy ', num2str(find(all(P2 < P2(:, i), 1))), '.']);
            P1_new(:, i) = [];
            P2_new(:, i) = [];
            % Dopo l'eliminazione di una colonna, ricomincia il loop
            % dalla prima colonna
            break;
        end
    end
end

function pureNashEquilibria(matrice1, matrice2)
    [num_righe1, num_colonne1] = size(matrice1);
    indici_minimi_matrice1 = [];

    % Trova gli indici minimi per le colonne della matrice1
    for j = 1:num_colonne1
        colonna_attuale = matrice1(:, j);
        valore_minimo = min(colonna_attuale);
        indici_minimi_colonna = find(colonna_attuale == valore_minimo);
        for i = 1:length(indici_minimi_colonna)
            indici_minimi_matrice1 = [indici_minimi_matrice1; indici_minimi_colonna(i), j];
        end
    end

    [num_righe2, num_colonne2] = size(matrice2);
    indici_minimi_matrice2 = [];

    % Trova gli indici minimi per le righe della matrice2
    for i = 1:num_righe2
        riga_attuale = matrice2(i, :);
        valore_minimo = min(riga_attuale);
        indici_minimi_riga = find(riga_attuale == valore_minimo);
        for j = 1:length(indici_minimi_riga)
            indici_minimi_matrice2 = [indici_minimi_matrice2; i, indici_minimi_riga(j)];
        end
    end

    % Calcola l'intersezione tra gli indici minimi delle due matrici
    intersezione = intersect(indici_minimi_matrice1, indici_minimi_matrice2, 'rows');

    if isempty(intersezione)
        disp('There are no pure Nash equilibria.');
    else
        disp('The pure Nash equilibria are:');
        disp(intersezione);
    end
end

function mixedNashEquilibria(matrice, var)
    syms x y; % Define x and y as symbolic variables

    mat_x = [x, 1 - x];
    mat_y = [y; 1 - y];

    risultato = mat_x * matrice * mat_y;
    disp(risultato)

    expr = simplify(expand(risultato));
    disp(expr)

    if var == 0
        collectedExpr = collect(expr, x);
        disp(collectedExpr)
    else
        collectedExpr = collect(expr, y);
        disp(collectedExpr)
    end
    
    analyzeExpression(expr, var)

end

function analyzeExpression(expr, var)
    syms x y;

    % Analyze based on variable
    if var == 0  % Analyze for x
        coeff = coeffs(expr, x);  % Get coefficient of x
        analyzeAndDisplay(coeff, 'x');  % Analyze and display solutions
    else  % Analyze for y
        coeff = coeffs(expr, y);  % Get coefficient of y
        analyzeAndDisplay(coeff, 'y');  % Analyze and display solutions
    end
end

function analyzeAndDisplay(coeff, var)
    term = sym(coeff(2));  % Convert to symbolic expression

    % Display based on the sign of the coefficient
    disp([char(term) ' > 0    ' var ' = 0']);
    disp([char(term) ' = 0    ' var ' E (0, 1)']);
    disp([char(term) ' < 0    ' var ' = 1']);

    disp(' ')
end

function mixedKKT(C1, C2)
    [m, n] = size(C1);
    
    H = [zeros(m, m), C1 + C2, ones(m, 1), zeros(m, 1); 
         (C1 + C2)', zeros(n, n), zeros(n, 1), ones(n, 1); 
         ones(1, m), zeros(1, n + 2); 
         zeros(1, m), ones(1, n), 0, 0];
     
    Ain = [-C2', zeros(n, n), zeros(n, 1), -ones(n, 1); 
           zeros(m, m), -C1, -ones(m, 1), zeros(m, 1)]; 
    bin = zeros(n + m, 1);
    
    Aeq = [ones(1, m), zeros(1, n + 2); 
           zeros(1, m), ones(1, n), 0, 0];

    numColonne = size(Aeq, 2);

    if mod(numColonne, 2) == 0
        primoValore = numColonne / 2;
        secondoValore = numColonne / 2;
    else
        primoValore = floor(numColonne / 2);
        secondoValore = ceil(numColonne / 2);
    end

    X0 = [rand(primoValore, 1); 10 - 20 * rand(secondoValore, 1)];
       
    beq = [1; 1];
    LB = [zeros(m + n, 1); -Inf; -Inf];
    UB = [ones(m + n, 1); Inf; Inf];
    
    [sol, fval, exitflag, output] = fmincon(@(X) 0.5 * X' * H * X, X0, Ain, bin, Aeq, beq, LB, UB);
    
    x = sol(1:m);
    disp("x:")
    disp(x)
    y = sol(m + 1:m + n);
    disp("y:")
    disp(y)
end
