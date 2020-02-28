close all;clc;clear all;
%%8by8 horizontal uniform%%
%%vertical SSH (modulated by vector a)%%
t1 = 1;t2 = 1;
r = -0.2;r0 = 0.2;k = 1;N = 8;a = [repmat([1.5 0.75],[1 N/2-1]) 1.2];b = 1*ones(1,N-1);
D = sparse(diag(-1i*r0 + repmat([0 1i*r],[1 N/2])) + diag(t1*b.*ones(1,N-1),1) + diag(t1*b.*ones(1,N-1),-1));
E = sparse(fliplr(diag(t2*ones(1,N))));F = sparse(N,N);
H = sparse([D E*a(1) F F F F F F;...
        E*a(1) D E*a(2) F F F F F;...
        F E*a(2) D E*a(3) F F F F;...
        F F E*a(3) D E*a(4) F F F;...
        F F F E*a(4) D E*a(5) F F;...
        F F F F E*a(5) D E*a(6) F;...
        F F F F F E*a(6) D E*a(7);...
        F F F F F F E*a(7) D]);
    
mov = VideoWriter('Movie0.avi');
mov.FrameRate = 2;
open(mov);
    %% middle 6 rows: S = [16 17 32 33 48 49];S = [12 21 28 37 44 53]; 
    %% all 8 rows: S = [1 16 17 32 33 48 49 64]; S = [5 12 21 28 37 44 53 60]; 
    %% first 7 rows: S = [1 16 17 32 33 48 49];
%%     S = [10 11 12 13 14 15 18 23 26 31 ]

    p = sparse(1,N^2);S = [8 10 11 21 22 28 29 35 ]; H = H + diag(p);D = 0.01*1i*eye(length(S));
%%
figure;

while true
    
    H(S,S)
    H(S,S) = H(S,S) + D;[A,V] = eig(full(H));Re = diag(real(V)/t1);Im = diag(imag(V)/t1);
    plot(Re,Im,'b*');hold on;lam = diag(V);xlim([-0.2 0.2]);xlabel('Real');ylabel('Imaginary');
    set(gca,'FontSize',14);
    set(gcf, 'Position', [00, 00, 400, 300]);
    if find(imag(lam)>0)
        break
    end
    k = k+1;
%     TEM = getframe;
%     writeVideo(mov,TEM.cdata);
    
end
close(mov);


Lasing = find(imag(diag(V))>0);Intensity = abs(A(:,Lasing)).*abs(A(:,Lasing));
Intensity = Intensity./max(Intensity);Phase = angle(A(:,Lasing));

figure;bar(Intensity);ylim([0 1]);set(gcf, 'Position', [00, 00, 400, 300]);

figure;bar(Phase);ylim([-pi pi]);set(gca,'ytick',[-pi,-pi/2,0,pi/2,pi]);
set(gca,'yticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'});set(gcf, 'Position', [00, 00, 400, 300]);

for k = 1:1:N
    if mod(k,2) RT(k,:) = Intensity((k-1)*N+(1:N));
    else
        RT(k,:) =  fliplr(Intensity((k-1)*N+(1:N))');
    end
end

for k = 1:1:N
    if mod(k,2)
        RA(k,:) = Phase((k-1)*N+(1:N));
    else
        RA(k,:) =  fliplr(Phase((k-1)*N+(1:N))');
    end
end


figure;imagesc(RT);colorbar;set(gcf, 'Position', [00, 00, 400, 300]);set(gca,'FontSize',14);

figure;imagesc(RA);colormap([0.5 0.5 0.5;0 1 1;0 1 0;1 0.5 0;0.5 0.5 0.5]);colorbar;   % X(:,2) = A(:,4);

