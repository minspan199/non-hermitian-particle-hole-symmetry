close all; clc; clear all;
t1 = 1; t2 = 1; r = -0.1; k = 1; N = 8; r0 = 0.5;
a = [repmat([0.574269153 1.595718388], [1 N/2-1]) 0.574269153]; b = ones(1, N - 1);
D = diag(-1i * r0 + repmat([0 1i*r], [1 N/2])) + diag(t1 * b .* ones(1, N - 1), 1) + diag(t1 * b .* ones(1, N - 1), -1);
E = fliplr(diag(t2 * ones(1, N)));
F = sparse(zeros(N));
H = sparse([D E * a(1) F F F F F F; ...
            E * a(1) D E * a(2) F F F F F; ...
            F E * a(2) D E * a(3) F F F F; ...
            F F E * a(3) D E * a(4) F F F; ...
            F F F E * a(4) D E * a(5) F F; ...
            F F F F E * a(5) D E * a(6) F; ...
            F F F F F E * a(6) D E * a(7); ...
            F F F F F F E * a(7) D; ]);
p = sparse(zeros(1, N^2));
s = zeros(1, N^2); D = 0.01 * 1i; S = 10; R = 316;
H = H + diag(p); H(S, S) = 0i;

figure;

while true
    H(S, S) = H(S, S) + D;
    [A, V] = eig(full(H));
    Re = diag(real(V) / t1);
    Im = diag(imag(V) / t1);
    plot(Re, Im, 'b*');
    hold on;
    xlabel('Real');
    ylabel('Imaginary');
    lam = diag(V);

    if find(imag(lam) > 0);
        break;
    end

    k = k + 1;
end

set(gcf, 'Position', [00, 00, 400, 300]);
xlim([-1 1]); set(gca, 'FontSize', 12);
Lasing = find(imag(diag(V)) > 0);
Intensity = abs(A(:, Lasing)) .* abs(A(:, Lasing));
Intensity = Intensity ./ max(Intensity);
Phase = angle(A(:, Lasing));

for k = 1:1:N

    if mod(k, 2)
        RT(k, :) = Intensity((k - 1) * N + (1:N));
    else
        RT(k, :) = fliplr(Intensity((k - 1) * N + (1:N))');
    end

end

for k = 1:1:N

    if mod(k, 2)
        RA(k, :) = Phase((k - 1) * N + (1:N));
    else
        RA(k, :) = fliplr(Phase((k - 1) * N + (1:N))');
    end

end

figure;
imagesc(RT);
colorbar;
set(gcf, 'Position', [00, 00, 400, 300]);
set(gca, 'fontsize', 12);
figure;
imagesc(RA);
colormap([0.5 0.5 0.5; 0 1 1; 0 1 0; 1 0.5 0; 0.5 0.5 0.5]);
colorbar;
set(gca, 'fontsize', 12);
set(gcf, 'Position', [00, 00, 400, 300]);
