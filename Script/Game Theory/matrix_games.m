function matrix_games()

    % 
    clc, clear all
    
    C = [
        5 4 3 5
        6 7 8 2
        6 3 4 4
    ];

    salva = C;
    stable = false;

    %Riduzione della Matrice
    while ~stable
        nuovaMatrice = strictly_dominated_strategies(C);
        stable = isequal(C, nuovaMatrice);

        C = nuovaMatrice;
    end

    C = salva;

    disp("Matrice Ridotta:");
    disp(nuovaMatrice)
    disp(" ")

    % Trova il minimo per ogni colonna (Player Righe)
    [num_righe_1, num_colonne_1] = size(C);
    result_matrix_1 = [];

    for j = 1:num_colonne_1
        colonna = C(:, j);
        minimo = min(colonna);
        indici_minimi = find(colonna == minimo);
        
        for i = 1:length(indici_minimi)
            result_matrix_1 = [result_matrix_1; indici_minimi(i), j];
        end
    end

    % Trova il massimo per ogni riga (Player Colonne)
    [num_righe_2, num_colonne_2] = size(C);
    result_matrix_2 = [];

    for i = 1:num_righe_2
        riga = C(i, :);
        massimo = max(riga);
        indici_massimi = find(riga == massimo);
        
        for j = 1:length(indici_massimi)
            result_matrix_2 = [result_matrix_2; i, indici_massimi(j)];
        end
    end

    % Confronto delle strategie di Nash
    intersezione = intersect(result_matrix_1, result_matrix_2, 'rows');

    if ~isempty(intersezione)
        disp('Strategie pure di Nash:');
        disp(intersezione);
    else
        disp('Non ci sono strategie pure di Nash');
        disp(' ');
    end

    disp("MIXED NASH EQUILIBRIUM")
    disp(" ")
    
    % Calculate m and n, the number of rows and columns of matrix C
    m = size(C, 1);
    n = size(C, 2);
    
    % Display the objective
    disp('min v')
    
    % Define the variable order (change the order as needed)
    variable_order = 1:m;
    
    % Dynamically create and display each inequality with the specified order
    for i = 1:n
        inequality = '';
        for j = variable_order
            % Append each term in the inequality to the string
            term = sprintf('%dx%d + ', C(j, i), j);
            inequality = [inequality term];
        end
        % Remove trailing plus and space, then add <= v
        inequality = [inequality(1:end-2) ' <= v'];
        disp(inequality)
    end
    
    % Display the equality constraint
    equality = '';
    for j = variable_order
        term = sprintf('x%d + ', j);
        equality = [equality term];
    end
    % Remove the trailing plus and space, then append '= 1'
    disp([equality(1:end-2) ' = 1'])
    
    % Display the non-negativity constraints
    nonnegativity = '';
    for j = variable_order
        nonnegativity = [nonnegativity sprintf('x%d, ', j)];
    end
    % Remove the trailing comma and space, then append '>= 0'
    disp([nonnegativity(1:end-2) ' >= 0'])

    if size(intersezione, 1) == 2 % Caso Particolare
        disp('Since the problem is linear then any convex combination alpha(x^1, v^1) + (1 - alpha)(x^2, v^2), alpha E [0, 1] is an optimal solution');
    else
        % Definisci c
        c = [zeros(m,1);1];
        
        % Definisci A
        A = [C', -ones(n,1)];
        b = zeros(n,1);
        
        % Definisci Aeq
        Aeq = [ones(1,m),0];
        beq = 1;
        
        % Definisci Ib
        Ib = [zeros(m,1);-inf];
        ub = [];
        
        % Calcola la soluzione
        [sol, Val,exitflag, output, lambda] = linprog(c, A,b, Aeq, beq, Ib, ub);
        disp('sol:')
        disp(sol(end, :)); % Stampa l'ultima riga di 'sol'
        disp('questa è l''ultima riga sol, che è v=w')
        
        % Estrai x
        x = sol(1:m);
        disp('x:')
        disp(x)
        
        % Estrai y
        y = lambda.ineqlin;
        disp('y:')
        disp(y)
    end
end

function [C_new] = strictly_dominated_strategies(C)
    % Inizializza le matrici P1_new e P2_new con le stesse dimensioni di P1 e P2
    C_new = C;

    % Elimina le righe in base alla condizione specificata
    for i = 1:size(C, 1)
        if any(all(C < C(i, :), 2))
            disp(['Strategy ', num2str(i), ' for Player 1 is strictly dominated by strategy ', num2str(find(all(C < C(i, :), 2))), '.']);
            C_new(i, :) = [];
            % Dopo l'eliminazione di una riga, ricomincia il loop
            % dalla prima riga
            break;
        end
    end

    % Elimina le colonne in base alla condizione specificata
    for i = 1:size(C, 2)
        if any(all(C > C(:, i), 1))
            disp(['Strategy ', num2str(i), ' for Player 2 is strictly dominated by strategy ', num2str(find(all(C > C(:, i), 1))), '.']);
            C_new(:, i) = [];
            % Dopo l'eliminazione di una colonna, ricomincia il loop
            % dalla prima colonna
            break;
        end
    end
end
