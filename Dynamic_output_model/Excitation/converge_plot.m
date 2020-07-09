%% Converge Plot
N_sample = [1500	2000	2500	3000	3500	4000	4500	5000];
standard = [0.618878578667 	0.699777049798 	0.282069927959 	0.362982274799 	0.018142010749 	0.018156235365 0.000195328600 0.000361520000];
S1  = [0.620744053	0.621397084	0.616614037	0.618204697	0.619667082	0.61890227	0.618706324	0.61871092];
S1t = [0.701345222	0.70134281	0.699903339	0.70127167	0.699830463	0.700301495	0.701091882	0.700577772];
S2  = [0.283109992	0.284265844	0.279754495	0.281395004	0.283354919	0.282488814	0.2812531	0.281783562];
S2t = [0.362921178	0.363125905	0.363138337	0.363804984	0.363680033	0.363564467	0.363405216	0.363669793];
S3  = [0.017790108	0.017891655	0.017709595	0.017775123	0.017692399	0.017997104	0.017978224	0.018220367];
S3t = [0.018477363	0.018408822	0.018414537	0.018445001	0.018396934	0.018497272	0.018479542	0.018485817];
SN = [0.000450349	0.000243655	0.000238625	0.000273508	0.00023817	0.000188279	0.000173835	0.000197112];
SNt= [0.000400212	0.000397265	0.000388144	0.000395763	0.000380847	0.000368852	0.000367208	0.00036174];

subplot(3,2,1)
plot(N_sample,S1);
hold on
plot(N_sample,S1,'r.','markersize',12)
plot(N_sample,ones(8,1)*standard(1));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_1$','Interpreter','latex')
axis([1500 5000 min(S1)*0.99 max(S1)*1.01 ]);

subplot(3,2,2)
plot(N_sample,S1t);
hold on
plot(N_sample,S1t,'r.','markersize',12)
plot(N_sample,ones(8,1)*standard(2));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_1^t$','Interpreter','latex')
axis([1500 5000 min(S1t)*0.99 max(S1t)*1.01 ]);

subplot(3,2,3)
plot(N_sample,S2);
hold on
plot(N_sample,S2,'r.','markersize',12)
plot(N_sample,ones(8,1)*standard(3));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_2$','Interpreter','latex')
axis([1500 5000 min(S2)*0.99 max(S2)*1.01 ]);

subplot(3,2,4)
plot(N_sample,S2t);
hold on
plot(N_sample,S2t,'r.','markersize',12)
plot(N_sample,ones(8,1)*standard(4));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_2^t$','Interpreter','latex')
axis([1500 5000 min(S2t)*0.99 max(S2t)*1.01 ]);

subplot(3,2,5)
plot(N_sample,S3);
hold on
plot(N_sample,S3,'r.','markersize',12)
plot(N_sample,ones(8,1)*standard(5));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_3$','Interpreter','latex')
axis([1500 5000 min(S3)*0.99 max(S3)*1.01 ]);

subplot(3,2,6)
plot(N_sample,S3t);
hold on
plot(N_sample,S3t,'r.','markersize',10)
plot(N_sample,ones(8,1)*standard(6));
xlabel('Sample Number','Interpreter','latex')
ylabel('Sensitivity Index','Interpreter','latex')
title('Converge of $S_3^t$','Interpreter','latex')
axis([1500 5000 min(S3t)*0.98 max(S3t)*1.02 ]);


