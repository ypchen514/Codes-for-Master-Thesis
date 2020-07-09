%% Converge Plot
clc
clear
close all
N_GSA = [1500 5000 10000 20000 30000 40000 60000 120000];

SS = [0.054978945	0.053507235	0.05337309	0.053638023	0.053866749	0.053896119	0.053845437	0.053868038;
0.651256349	0.640955191	0.640898561	0.640683392	0.640953722	0.640931972	0.640877626	0.640788606;
0.202502633	0.199752902	0.197952102	0.198816813	0.198542293	0.198743226	0.1983113	0.197821085;
0.001502589	0.001248351	0.000968383	0.000839445	0.000707612	0.000822391	0.000754226	0.000574512;
0.073129192	0.071254071	0.071497824	0.071457282	0.071734854	0.071790035	0.071771788	0.071774691;
0.001392922	0.000969563	0.00058079	0.000218847	0.000292812	0.000407426	0.00034091	0.000329209;
0.064869981	0.065036078	0.06500274	0.065033987	0.065082392	0.065171443	0.065174969	0.065141832;
0.665055791	0.668306524	0.668308739	0.669245861	0.669391298	0.6697163	0.669360549	0.669475604;
0.226953724	0.221634005	0.221013731	0.220797631	0.221334873	0.22150044	0.221379882	0.221039055;
0.003432611	0.003348177	0.00332325	0.003311564	0.003320788	0.003316834	0.003347193	0.003325993;
0.077877609	0.078147449	0.078412207	0.078369042	0.07838644	0.078488706	0.078447372	0.078509068;
0.003058564	0.003093951	0.003031336	0.003015342	0.002986154	0.003000161	0.002997673	0.00299163];


% Dlane fusion data
S1 = SS(1,:);
S2 = SS(2,:);
S3 = SS(3,:);
S4 = SS(4,:);
S5 = SS(5,:);
S6 = SS(6,:);
St1 = SS(7,:);
St2 = SS(8,:);
St3 = SS(9,:);
St4 = SS(10,:);
St5 = SS(11,:);
St6 = SS(12,:);

figure
%1m
subplot(6,2,1)
plot(log10(N_GSA),S1);
hold on
plot(log10(N_GSA),S1,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S1)-0.001) (max(S1)+0.001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_1$','Interpreter','latex');
%1t
subplot(6,2,2)
plot(log10(N_GSA),St1);
hold on
plot(log10(N_GSA),St1,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St1)-0.001) (max(St1)+0.001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_1$','Interpreter','latex');
%2m
subplot(6,2,3)
plot(log10(N_GSA),S2);
hold on
plot(log10(N_GSA),S2,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S2)-0.005) (max(S2)+0.001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_2$','Interpreter','latex');
%2t
subplot(6,2,4)
plot(log10(N_GSA),St2);
hold on
plot(log10(N_GSA),St2,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St2)-0.001) (max(St2)+0.002)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_2$','Interpreter','latex');
%3m
subplot(6,2,5)
plot(log10(N_GSA),S3);
hold on
plot(log10(N_GSA),S3,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S3)-0.005) (max(S3)+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_3$','Interpreter','latex');
%3t
subplot(6,2,6)
plot(log10(N_GSA),St3);
hold on
plot(log10(N_GSA),St3,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St3)-0.005) (max(St3)+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_3$','Interpreter','latex');
%4m
subplot(6,2,7)
plot(log10(N_GSA),S4);
hold on
plot(log10(N_GSA),S4,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S4)-0.005) (max(S4)+0.0001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_4$','Interpreter','latex');
%4t
subplot(6,2,8)
plot(log10(N_GSA),St4);
hold on
plot(log10(N_GSA),St4,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St4)-0.0002) (max(St4)+0.0001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_4^t$','Interpreter','latex');

%5m
subplot(6,2,9)
plot(log10(N_GSA),S5);
hold on
plot(log10(N_GSA),S5,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S5)-0.001) (max(S5)+0.001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_5$','Interpreter','latex');
%5t
subplot(6,2,10)
plot(log10(N_GSA),St5);
hold on
plot(log10(N_GSA),St5,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St5)-0.001) (max(St5)+0.001)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_5$','Interpreter','latex');

%6m
subplot(6,2,11)
plot(log10(N_GSA),S6);
hold on
plot(log10(N_GSA),S6,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(S6)-0.005) (max(S6)+0.0005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_6$','Interpreter','latex');
%6t
subplot(6,2,12)
plot(log10(N_GSA),St6);
hold on
plot(log10(N_GSA),St6,'r*');
axis([log10(N_GSA(1)) log10(N_GSA(8)) max(0,min(St6)-0.0002) (max(St6)+0.0002)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_6$','Interpreter','latex');






