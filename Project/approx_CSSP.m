function R = approx_CSSP(A, k, ITER_NUM)
% Input Arguments:
% A: the matrix
% k: the number of columns to choose
% ITER_NUM: the number of repetitive iterations
%
% Output Arguments:
% R: the set of column indices
%
% Reference:
% Bruno Ordozgoiti, Sandra Gomez Canaval, Alberto Mozo, "A Fast Iterative 
% Algorithm for Improved Unsupervised Feature Selection". IEEE International 
% Conference on Data Mining, 2016.

    if (nargin < 3)
        ITER_NUM = Inf;
    end
    [m, n] = size(A);    
    % Dealing with huge matrices, Theorem 2
    if(m > n)
        [U,S,V] = svd(A, 'econ');
        A = S * V';
    end
    subset = randsample(n, k, false);
    former_subset = subset;
    former_norm = 0;
    C = A(:, subset);
    pinv_C = pinv(C);
    E = A - C*pinv_C*A;
    
    % Dealing with huge matrices, Theorem 3
    if(n > m)
        [U,S,V] = svd(E, 'econ');
        F = S * S * V';
        G = S * V';
    else
        F = E' * E;
        G = E;
    end
    f = sum(F.^2);
    g = sum(G.^2);    
    
    % substitution
    flag = 1;
    while flag
        for i = 1:k
            j = subset(i);
            
            % Proposition 1
            temp_C = C;
            % Zero out column i
            temp_C(:, i) = 0;
            rou = (pinv_C(i, :))';
            pinv_temp_C = pinv_C - 1/(rou'*rou)*pinv_C*(rou*rou');
            
            % Proposition 2
            S1 = C(:, i) * rou' * A;
            S2 = 1/(rou'*rou) * temp_C * pinv_C * (rou*rou') * A;
            temp_E = E + S1 + S2;
            
            % Proposition 3
            delta = temp_E(:, j)' * temp_E;
            gamma = E' * E * delta';
            temp_f = f + norm(delta)^2*(delta.*delta)/(delta(j)^2) + 2*(gamma'.*delta)/delta(j);
            temp_g = g + (delta.*delta)/delta(j);
            
            % Theorem 1
            scores = temp_f./ temp_g;
            for h = subset'
                if h ~= subset(i)
                     scores(h) = 0;
                end
            end
            scores(temp_g < 1e-5) = 0;
            
            % find w
            [~, w] = max(scores);
            delta = temp_E(:, w)' * temp_E;
            gamma = temp_E' * temp_E * delta';
            f = temp_f + norm(delta)^2*(delta.*delta)/(delta(w)^2) - 2*(gamma'.*delta)/delta(w);
            g = temp_g - (delta.*delta)/delta(w);
            
            % Proposition 4
            omega = temp_E(:, w);
            temp = 1 / (norm(omega)^2) * (pinv_temp_C*A(:, w)*omega');
            temp(i, :) = -1 / (norm(omega)^2) * omega';
            pinv_C = pinv_temp_C - temp;
            
            % update
            E = temp_E - omega*omega'*temp_E/(norm(omega)^2);
            subset(i) = w;
            C = A(:,subset);
        end
        % stop criteria
        if(ITER_NUM > 1)
            if(norm(C*pinv(C)*A, 'fro') > former_norm)
                former_norm = norm(C*pinv(C)*A, 'fro');
                former_subset = subset;
            else
                subset = former_subset;
                flag = 0;
            end
        else
            flag = 0;
        end
    end
    % output R
    R = subset;
end