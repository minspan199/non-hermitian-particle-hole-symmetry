close all;
clc;
clear all;
t1 = 1; t2 = 1;
r = -0.2; k = 1;
N = 6; r0 = 0.5;
a = [repmat([0.596870221 1.592352575], [1 N/2-1]) 1.592352575];
b = ones(1, N - 1) .* [repmat([0.596870221 1.592352575], [1 N/2-1]) 1.592352575];
D = diag(-1i * r0 + repmat([0 1i*r], [1 N/2])) + diag(t1 * b .* ones(1, N - 1), 1) + diag(t1 * b .* ones(1, N - 1), -1);
E = fliplr(diag(t2 * ones(1, N)));
F = sparse(zeros(N));

H = sparse([D E * a(1) F F F F; ...
            E * a(1) D E * a(2) F F F; ...
            F E * a(2) D E * a(3) F F; ...
            F F E * a(3) D E * a(4) F; ...
            F F F E * a(4) D E * a(5); ...
            F F F F E * a(5) D]);
%% open circle 3by3: S = [1 12 13 14 15 10 3];
%% backward diagonal: S = [1 12 11 14 15 22 21 28 29 32 31];
%% first column: S = [1 12 13 24 25 36];
%% S = [10 11 12 13 14 15 18 23 26 31];
%% S = [11 14 15 22 21 28 29];
%% S = [11 14 15 22 21]
% p = sparse(1,N^2);S = [8 10 11 21 22 28 29 35 ]; H = H + diag(p);D = 0.01*1i*eye(length(S));

% S = [1 12 13 24 25 36];

Mt(36) = 0;

for i = 1:1:36
    i

    Ht = H;
    D = 0.01 * 1i * eye(length(i));

    while true
        Ht(i, i) = Ht(i, i) + D;
        [A, V] = eig(full(Ht));
        Re = diag(real(V) / t1);
        Im = diag(imag(V) / t1);
        lam = diag(V);

        if find(imag(lam) > 0)
            break
        end

        k = k + 1;
    end

    Lasing = find(imag(diag(V)) > 0);

    if (length(Lasing)==1&&real(lam(Lasing)) < 1e3)
         real(lam(Lasing))
        Mt(i) = 1;
    else
        Mt(i) = 0;
    end

end


St = [[1,0,0,1,0,1];[0,1,1,1,1,1];[1,1,1,1,1,0];[0,1,1,1,1,1];[1,1,1,1,1,0];[1,1,1,1,1,1]];
figure;
imagesc(St);

