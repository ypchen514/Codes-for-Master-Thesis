%% Direct Optimization Algorithm - - Main
% main file
clc
clear

tic
% setting
options.conTol = 1e-3;
options.display = 0;
fileInfo.fName = 'obj';
% algorithm
result = UMDIRECT(fileInfo,[5,0,0],[10,0.5,10], options);
% result
x_best = result.xBest;
y_best = result.fBest;
toc

history = result.rect.x;

parfor i = 1:length(history)
    [S_overall,var] = F_sensitivity(history(i,:));
    S_1(i) = S_overall(1,1);
    St_1(i) = S_overall(1,2);
    S_2(i) = S_overall(2,1);
    St_2(i) = S_overall(2,2);
    S_3(i) = S_overall(3,1);
    St_3(i) = S_overall(3,2);
    Var(i) = var; 
end

subplot(3,3,1)
plot(0.6480,0.7254,'o','markersize',14)
plot(0.6480,0.7254,'o','markersize',14)
plot(S_1,St_1,'g.','markersize',12)
xlabel('$S_1$','interpreter','latex','Fontsize',12);
ylabel('$S^t_1$','interpreter','latex','Fontsize',12);
hold on
plot(0.6480,0.7254,'r.','markersize',14)
title('$S_1$ vs $S_1^t$','interpreter','latex','Fontsize',16);
axis([0.4,0.7,0,1])

subplot(3,3,2)
plot(0.6469,0.2756+0.00041538,'r.','markersize',14)
plot(S_1,S_2+S_3,'g.','markersize',12)
hold on
plot(0.6469,0.2756+0.00041538,'r.','markersize',14)
xlabel('$S_1$','interpreter','latex','Fontsize',12);
ylabel('$S_2+S_3$','interpreter','latex','Fontsize',12);
title('$S_1$ vs $(S_2+S_3)$','interpreter','latex','Fontsize',16);

subplot(3,3,3)
plot(0.6469,11977,'r.','markersize',14)
plot(S_1,Var,'g.','markersize',12)
hold on
plot(0.6469,11977,'r.','markersize',14)
xlabel('$S_1$','interpreter','latex','Fontsize',12);
ylabel('$Var$','interpreter','latex','Fontsize',12);
title('$S_1$ vs $Var$','interpreter','latex','Fontsize',16);
% 
% subplot(3,3,4)
% plot(S_1,St_1,'g.','markersize',12)
% xlabel('$S_1$','interpreter','latex','Fontsize',12);
% ylabel('$S^t_1$','interpreter','latex','Fontsize',12);
% hold on
% plot(0.6480,0.7254,'r.','markersize',14)
% title('$S_1$ vs $S_1^t$','interpreter','latex','Fontsize',16);
% 
% subplot(3,3,2)
% plot(S_1,S_2+S_3,'g.','markersize',12)
% hold on
% plot(0.6469,0.2756+0.00041538,'r.','markersize',14)
% xlabel('$S_1$','interpreter','latex','Fontsize',12);
% ylabel('$S_2+S_3$','interpreter','latex','Fontsize',12);
% title('$S_1$ vs $S_2+S_3$','interpreter','latex','Fontsize',16);
% 
% subplot(3,3,3)
% plot(S_1,Var,'g.','markersize',12)
% hold on
% plot(0.6469,11977,'r.','markersize',14)
% xlabel('$S_1$','interpreter','latex','Fontsize',12);
% ylabel('$Var$','interpreter','latex','Fontsize',12);
% title('$S_1$ vs $Var$','interpreter','latex','Fontsize',16);

% % 
% subplot(3,3,4)
% plot(S_2,St_2,'g.','markersize',12)
% xlabel('$S_2$','interpreter','latex','Fontsize',12);
% ylabel('$S^t_2$','interpreter','latex','Fontsize',12);
% hold on
% plot(0.4809,0.5271,'r.','markersize',14)
% title('$S_2$ vs $S_2^t$','interpreter','latex','Fontsize',16);
% % 
% subplot(3,3,5)
% plot(S_2,S_1+S_3,'g.','markersize',12)
% hold on
% plot(0.4809,0.4329+0.0425,'r.','markersize',14)
% xlabel('$S_2$','interpreter','latex','Fontsize',12);
% ylabel('$S_2+S_3$','interpreter','latex','Fontsize',12);
% title('$S_2$ vs $(S_1+S_3)$','interpreter','latex','Fontsize',16);
% % 
% subplot(3,3,6)
% plot(S_2,Var,'g.','markersize',12)
% hold on
% plot(0.4809,9661,'r.','markersize',14)
% xlabel('$S_2$','interpreter','latex','Fontsize',12);
% ylabel('$Var$','interpreter','latex','Fontsize',12);
% title('$S_2$ vs $Var$','interpreter','latex','Fontsize',16);
% % 
subplot(3,3,7)
hold on
plot(0.1065,0.1067,'bo','markersize',6)
plot(0.1678,0.1678,'r.','markersize',14)
plot(S_3,St_3,'g.','markersize',12)
xlabel('$S_3$','interpreter','latex','Fontsize',12);
ylabel('$S^t_3$','interpreter','latex','Fontsize',12);
plot(0.1678,0.1678,'r.','markersize',14)
plot(0.1065,0.1067,'bo','markersize',8)
legend('Kriging','Direct','interpreter','latex')
title('$S_3$ vs $S_3^t$','interpreter','latex','Fontsize',16);
axis([0,0.5,0,0.5])
% 
subplot(3,3,8)
hold on
plot(0.1067,0.4969+0.3414,'bo','markersize',6)
plot(0.1678,0.3418+0.4461,'r.','markersize',14)
plot(S_3,S_1+S_2,'g.','markersize',12)
plot(0.1678,0.3418+0.4461,'r.','markersize',14)
plot(0.1067,0.4969+0.3414,'bo','markersize',6)
legend('Kriging','Direct','interpreter','latex')
xlabel('$S_3$','interpreter','latex','Fontsize',12);
ylabel('$S_1+S_2$','interpreter','latex','Fontsize',12);
title('$S_3$ vs $(S_1+S_2)$','interpreter','latex','Fontsize',16);
axis([0,0.2,0.6,1])
% 
subplot(3,3,9)
hold off

plot(0.167,4247.5,'bo','markersize',6)
hold on
plot(0.1678,3454,'r.','markersize',14)
plot(S_3,Var,'g.','markersize',12)
plot(0.1678,3454,'r.','markersize',14)
plot(0.1067,42475,'bo','markersize',6)
xlabel('$S_3$','interpreter','latex','Fontsize',12);
ylabel('$Var$','interpreter','latex','Fontsize',12);
title('$S_3$ vs $Var$','interpreter','latex','Fontsize',16);
legend('Kriging','Direct','interpreter','latex')
