% x = [0.05288 	0.05131 	0.05496 	0.05163 	0.05259 	0.05229 	0.05293 
% 0.22710 	0.21636 	0.21393 	0.21608 	0.21669 	0.21704 	0.21681 
% 0.07060 	0.07424 	0.07441 	0.07216 	0.07297 	0.07291 	0.07406 
% 0.06719 	0.06832 	0.06970 	0.06922 	0.06901 	0.06897 	0.06920 
% 0.05750 	0.05335 	0.05608 	0.05666 	0.05573 	0.05593 	0.05581 
% 0.04496 	0.05106 	0.05332 	0.05228 	0.05204 	0.05250 	0.05267 
% ];
% 
% y = [0.09007 	0.08789 	0.09066 	0.08871 	0.08943 	0.08928 	0.08979 
% 0.17785 	0.17031 	0.17050 	0.17222 	0.17304 	0.17314 	0.17346 
% 0.07877 	0.08438 	0.08405 	0.08300 	0.08351 	0.08355 	0.08412 
% 0.11450 	0.11336 	0.11503 	0.11429 	0.11418 	0.11424 	0.11445 
% 0.09985 	0.09397 	0.09519 	0.09540 	0.09462 	0.09437 	0.09424 
% 0.08684 	0.08893 	0.09076 	0.08987 	0.08963 	0.09009 	0.09018 
% ];
% 
% yaw = [0.07087 	0.06935 	0.07299 	0.06970 	0.07075 	0.07044 	0.07105 
% 0.21551 	0.20591 	0.20420 	0.20624 	0.20692 	0.20721 	0.20719 
% 0.07637 	0.08002 	0.08025 	0.07813 	0.07894 	0.07882 	0.07979 
% 0.07239 	0.07337 	0.07490 	0.07418 	0.07393 	0.07391 	0.07412 
% 0.05719 	0.05338 	0.05560 	0.05636 	0.05535 	0.05557 	0.05572 
% 0.05094 	0.05650 	0.05853 	0.05764 	0.05739 	0.05781 	0.05801 
% ];
% 
fusion = [0.06207 	0.06035 	0.06379 	0.06079 	0.06169 	0.06143 	0.06204 
0.21496 	0.20501 	0.20322 	0.20527 	0.20593 	0.20622 	0.20613 
0.07262 	0.07675 	0.07680 	0.07484 	0.07558 	0.07554 	0.07655 
0.07884 	0.07941 	0.08087 	0.08032 	0.08014 	0.08012 	0.08034 
0.06792 	0.06334 	0.06570 	0.06619 	0.06530 	0.06538 	0.06527 
0.05528 	0.06039 	0.06254 	0.06154 	0.06130 	0.06176 	0.06191 
];
% 
% 
% x_t = [0.19795 	0.19922 	0.19814 	0.19849 	0.19829 	0.19797 	0.19776 
% 0.52047 	0.53338 	0.53549 	0.53441 	0.53455 	0.53427 	0.53377 
% 0.23971 	0.24006 	0.24044 	0.24013 	0.23806 	0.23789 	0.23802 
% 0.21281 	0.21240 	0.21439 	0.21738 	0.21652 	0.21702 	0.21678 
% 0.21670 	0.22019 	0.21983 	0.21956 	0.21941 	0.21965 	0.21943 
% 0.21000 	0.21251 	0.21195 	0.21301 	0.21313 	0.21294 	0.21276 
% ];
% 
% y_t = [0.19851 	0.19986 	0.19911 	0.19893 	0.19898 	0.19874 	0.19865 
% 0.35806 	0.36446 	0.36655 	0.36553 	0.36563 	0.36550 	0.36524 
% 0.19609 	0.19463 	0.19409 	0.19361 	0.19268 	0.19264 	0.19269 
% 0.23479 	0.23508 	0.23549 	0.23703 	0.23661 	0.23681 	0.23677 
% 0.22187 	0.22143 	0.22090 	0.22033 	0.22020 	0.22036 	0.22020 
% 0.21225 	0.21221 	0.21174 	0.21224 	0.21255 	0.21215 	0.21209 
% ];
% 
% yaw_t = [0.21661 	0.21778 	0.21685 	0.21708 	0.21689 	0.21662 	0.21640 
% 0.48499 	0.49712 	0.49897 	0.49806 	0.49832 	0.49803 	0.49757 
% 0.23465 	0.23500 	0.23547 	0.23509 	0.23329 	0.23315 	0.23326 
% 0.21124 	0.21091 	0.21252 	0.21525 	0.21446 	0.21491 	0.21471 
% 0.20805 	0.21052 	0.21054 	0.21032 	0.21019 	0.21037 	0.21016 
% 0.20494 	0.20779 	0.20716 	0.20805 	0.20822 	0.20802 	0.20787 
% ];

fusion_t = [0.19813 	0.19942 	0.19842 	0.19864 	0.19850 	0.19820 	0.19802 
0.48044 	0.49175 	0.49384 	0.49278 	0.49291 	0.49267 	0.49223 
0.22897 	0.22887 	0.22903 	0.22867 	0.22689 	0.22675 	0.22686 
0.21821 	0.21798 	0.21958 	0.22221 	0.22146 	0.22188 	0.22170 
0.21795 	0.22047 	0.22008 	0.21973 	0.21958 	0.21981 	0.21960 
0.21055 	0.21243 	0.21188 	0.21281 	0.21297 	0.21273 	0.21258 
];

N_sample = [5000 	10000 	30000 	60000 	100000 	120000 	150000 
];

%% plot
% figure
% %1m
% subplot(6,2,1)
% plot(log10(N_sample),x(1,:));
% hold on
% plot(log10(N_sample),x(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,mean(x(1,:))-0.005) (mean(x(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_1$','Interpreter','latex');
% %1t
% subplot(6,2,2)
% plot(log10(N_sample),x_t(1,:));
% hold on
% plot(log10(N_sample),x_t(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(1,:))-0.005) (max(x_t(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_1$','Interpreter','latex');
% %2m
% subplot(6,2,3)
% plot(log10(N_sample),x(2,:));
% hold on
% plot(log10(N_sample),x(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x(2,:))-0.005) (max(x(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_2$','Interpreter','latex');
% %2t
% subplot(6,2,4)
% plot(log10(N_sample),x_t(2,:));
% hold on
% plot(log10(N_sample),x_t(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(2,:))-0.005) (max(x_t(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_2$','Interpreter','latex');
% %3m
% subplot(6,2,5)
% plot(log10(N_sample),x(3,:));
% hold on
% plot(log10(N_sample),x(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x(3,:))-0.005) (max(x(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_3$','Interpreter','latex');
% %3t
% subplot(6,2,6)
% plot(log10(N_sample),x_t(3,:));
% hold on
% plot(log10(N_sample),x_t(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(3,:))-0.005) (max(x_t(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_3$','Interpreter','latex');
% %4m
% subplot(6,2,7)
% plot(log10(N_sample),x(4,:));
% hold on
% plot(log10(N_sample),x(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x(4,:))-0.005) (max(x(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_4$','Interpreter','latex');
% %4t
% subplot(6,2,8)
% plot(log10(N_sample),x_t(4,:));
% hold on
% plot(log10(N_sample),x_t(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(4,:))-0.005) (max(x_t(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% 
% %5m
% subplot(6,2,9)
% plot(log10(N_sample),x(5,:));
% hold on
% plot(log10(N_sample),x(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x(5,:))-0.005) (max(x(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_5$','Interpreter','latex');
% %5t
% subplot(6,2,10)
% plot(log10(N_sample),x_t(5,:));
% hold on
% plot(log10(N_sample),x_t(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(5,:))-0.005) (max(x_t(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_5$','Interpreter','latex');
% 
% %6m
% subplot(6,2,11)
% plot(log10(N_sample),x(6,:));
% hold on
% plot(log10(N_sample),x(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x(6,:))-0.005) (max(x(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_6$','Interpreter','latex');
% %6t
% subplot(6,2,12)
% plot(log10(N_sample),x_t(6,:));
% hold on
% plot(log10(N_sample),x_t(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(x_t(6,:))-0.005) (max(x_t(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_6$','Interpreter','latex');
% 
% figure
% %1m
% subplot(6,2,1)
% plot(log10(N_sample),y(1,:));
% hold on
% plot(log10(N_sample),y(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,mean(y(1,:))-0.005) (mean(y(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_1$','Interpreter','latex');
% %1t
% subplot(6,2,2)
% plot(log10(N_sample),y_t(1,:));
% hold on
% plot(log10(N_sample),y_t(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(1,:))-0.005) (max(y_t(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_1$','Interpreter','latex');
% %2m
% subplot(6,2,3)
% plot(log10(N_sample),y(2,:));
% hold on
% plot(log10(N_sample),y(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y(2,:))-0.005) (max(y(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_2$','Interpreter','latex');
% %2t
% subplot(6,2,4)
% plot(log10(N_sample),y_t(2,:));
% hold on
% plot(log10(N_sample),y_t(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(2,:))-0.005) (max(y_t(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_2$','Interpreter','latex');
% %3m
% subplot(6,2,5)
% plot(log10(N_sample),y(3,:));
% hold on
% plot(log10(N_sample),y(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y(3,:))-0.005) (max(y(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_3$','Interpreter','latex');
% %3t
% subplot(6,2,6)
% plot(log10(N_sample),y_t(3,:));
% hold on
% plot(log10(N_sample),y_t(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(3,:))-0.005) (max(y_t(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_3$','Interpreter','latex');
% %4m
% subplot(6,2,7)
% plot(log10(N_sample),y(4,:));
% hold on
% plot(log10(N_sample),y(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y(4,:))-0.005) (max(y(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_4$','Interpreter','latex');
% %4t
% subplot(6,2,8)
% plot(log10(N_sample),y_t(4,:));
% hold on
% plot(log10(N_sample),y_t(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(4,:))-0.005) (max(y_t(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% 
% %5m
% subplot(6,2,9)
% plot(log10(N_sample),y(5,:));
% hold on
% plot(log10(N_sample),y(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y(5,:))-0.005) (max(y(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_5$','Interpreter','latex');
% %5t
% subplot(6,2,10)
% plot(log10(N_sample),y_t(5,:));
% hold on
% plot(log10(N_sample),y_t(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(5,:))-0.005) (max(y_t(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_5$','Interpreter','latex');
% 
% %6m
% subplot(6,2,11)
% plot(log10(N_sample),y(6,:));
% hold on
% plot(log10(N_sample),y(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y(6,:))-0.005) (max(y(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_6$','Interpreter','latex');
% %6t
% subplot(6,2,12)
% plot(log10(N_sample),y_t(6,:));
% hold on
% plot(log10(N_sample),y_t(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(y_t(6,:))-0.005) (max(y_t(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_6$','Interpreter','latex');
% 
% figure
% %1m
% subplot(6,2,1)
% plot(log10(N_sample),yaw(1,:));
% hold on
% plot(log10(N_sample),yaw(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,mean(yaw(1,:))-0.005) (mean(yaw(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_1$','Interpreter','latex');
% %1t
% subplot(6,2,2)
% plot(log10(N_sample),yaw_t(1,:));
% hold on
% plot(log10(N_sample),yaw_t(1,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(1,:))-0.005) (max(yaw_t(1,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_1$','Interpreter','latex');
% %2m
% subplot(6,2,3)
% plot(log10(N_sample),yaw(2,:));
% hold on
% plot(log10(N_sample),yaw(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw(2,:))-0.005) (max(yaw(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_2$','Interpreter','latex');
% %2t
% subplot(6,2,4)
% plot(log10(N_sample),yaw_t(2,:));
% hold on
% plot(log10(N_sample),yaw_t(2,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(2,:))-0.005) (max(yaw_t(2,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_2$','Interpreter','latex');
% %3m
% subplot(6,2,5)
% plot(log10(N_sample),yaw(3,:));
% hold on
% plot(log10(N_sample),yaw(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw(3,:))-0.005) (max(yaw(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_3$','Interpreter','latex');
% %3t
% subplot(6,2,6)
% plot(log10(N_sample),yaw_t(3,:));
% hold on
% plot(log10(N_sample),yaw_t(3,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(3,:))-0.005) (max(yaw_t(3,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_3$','Interpreter','latex');
% %4m
% subplot(6,2,7)
% plot(log10(N_sample),yaw(4,:));
% hold on
% plot(log10(N_sample),yaw(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw(4,:))-0.005) (max(yaw(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_4$','Interpreter','latex');
% %4t
% subplot(6,2,8)
% plot(log10(N_sample),yaw_t(4,:));
% hold on
% plot(log10(N_sample),yaw_t(4,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(4,:))-0.005) (max(yaw_t(4,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% 
% %5m
% subplot(6,2,9)
% plot(log10(N_sample),yaw(5,:));
% hold on
% plot(log10(N_sample),yaw(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw(5,:))-0.005) (max(yaw(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_5$','Interpreter','latex');
% %5t
% subplot(6,2,10)
% plot(log10(N_sample),yaw_t(5,:));
% hold on
% plot(log10(N_sample),yaw_t(5,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(5,:))-0.005) (max(yaw_t(5,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_5$','Interpreter','latex');
% 
% %6m
% subplot(6,2,11)
% plot(log10(N_sample),yaw(6,:));
% hold on
% plot(log10(N_sample),yaw(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw(6,:))-0.005) (max(yaw(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S_6$','Interpreter','latex');
% %6t
% subplot(6,2,12)
% plot(log10(N_sample),yaw_t(6,:));
% hold on
% plot(log10(N_sample),yaw_t(6,:),'r*');
% axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(yaw_t(6,:))-0.005) (max(yaw_t(6,:))+0.005)])
% xlabel('$\log (N_{GSA})$','Interpreter','latex');
% ylabel('Sensitivity','Interpreter','latex');
% title('Covergence of $S^t_6$','Interpreter','latex');

figure
%1m
subplot(6,2,1)
plot(log10(N_sample),fusion(1,:));
hold on
plot(log10(N_sample),fusion(1,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,mean(fusion(1,:))-0.005) (mean(fusion(1,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_1$','Interpreter','latex');
%1t
subplot(6,2,2)
plot(log10(N_sample),fusion_t(1,:));
hold on
plot(log10(N_sample),fusion_t(1,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(1,:))-0.005) (max(fusion_t(1,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_1$','Interpreter','latex');
%2m
subplot(6,2,3)
plot(log10(N_sample),fusion(2,:));
hold on
plot(log10(N_sample),fusion(2,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion(2,:))-0.005) (max(fusion(2,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_2$','Interpreter','latex');
%2t
subplot(6,2,4)
plot(log10(N_sample),fusion_t(2,:));
hold on
plot(log10(N_sample),fusion_t(2,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(2,:))-0.005) (max(fusion_t(2,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_2$','Interpreter','latex');
%3m
subplot(6,2,5)
plot(log10(N_sample),fusion(3,:));
hold on
plot(log10(N_sample),fusion(3,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion(3,:))-0.005) (max(fusion(3,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_3$','Interpreter','latex');
%3t
subplot(6,2,6)
plot(log10(N_sample),fusion_t(3,:));
hold on
plot(log10(N_sample),fusion_t(3,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(3,:))-0.005) (max(fusion_t(3,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_3$','Interpreter','latex');
%4m
subplot(6,2,7)
plot(log10(N_sample),fusion(4,:));
hold on
plot(log10(N_sample),fusion(4,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion(4,:))-0.005) (max(fusion(4,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_4$','Interpreter','latex');
%4t
subplot(6,2,8)
plot(log10(N_sample),fusion_t(4,:));
hold on
plot(log10(N_sample),fusion_t(4,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(4,:))-0.005) (max(fusion_t(4,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');

%5m
subplot(6,2,9)
plot(log10(N_sample),fusion(5,:));
hold on
plot(log10(N_sample),fusion(5,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion(5,:))-0.005) (max(fusion(5,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_5$','Interpreter','latex');
%5t
subplot(6,2,10)
plot(log10(N_sample),fusion_t(5,:));
hold on
plot(log10(N_sample),fusion_t(5,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(5,:))-0.005) (max(fusion_t(5,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_5$','Interpreter','latex');

%6m
subplot(6,2,11)
plot(log10(N_sample),fusion(6,:));
hold on
plot(log10(N_sample),fusion(6,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion(6,:))-0.005) (max(fusion(6,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S_6$','Interpreter','latex');
%6t
subplot(6,2,12)
plot(log10(N_sample),fusion_t(6,:));
hold on
plot(log10(N_sample),fusion_t(6,:),'r*');
axis([log10(N_sample(1)) log10(N_sample(7)) max(0,min(fusion_t(6,:))-0.005) (max(fusion_t(6,:))+0.005)])
xlabel('$\log (N_{GSA})$','Interpreter','latex');
ylabel('Sensitivity','Interpreter','latex');
title('Covergence of $S^t_6$','Interpreter','latex');

