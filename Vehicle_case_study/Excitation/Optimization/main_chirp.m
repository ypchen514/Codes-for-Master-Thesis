%% Kriging Optimization--chirp
clear
clc


global S1_krig S2_krig S3_krig S4_krig S5_krig S6_krig St1_krig St2_krig St3_krig St4_krig St5_krig St6_krig Var_krig count
global S1_est S2_est S3_est S4_est S5_est S6_est St1_est St2_est St3_est St4_est St5_est St6_est var_est

count = 0;


input = [0.05	0.05	0.05	0.05	0.05	0.05	0.1	0.1	0.1	0.1	0.1	0.1	0.15	0.15	0.15	0.15	0.15	0.15	0.2	0.2	0.2	0.2	0.2	0.2	0.25	0.25	0.25	0.25	0.25	0.25;
0.5	0.6	0.7	0.8	0.9	1	0.5	0.6	0.7	0.8	0.9	1	0.5	0.6	0.7	0.8	0.9	1	0.5	0.6	0.7	0.8	0.9	1	0.5	0.6	0.7	0.8	0.9	1]';
lb = [0.05,0.5];
ub = [0.25,1];
fittype = 1;
SCFtype = 1;
max_variance = 1;
% S1
output_S1 = [0.037888044	0.040348768	0.041971832	0.042925081	0.043519041	0.050480865	0.070614227	0.067482423	0.061552663	0.053743411	0.045276182	0.037112691	0.064009279	0.052323118	0.041106251	0.032572527	0.029895691	0.036771451	0.056221063	0.041427541	0.034073938	0.042708322	0.069448057	0.111800378	0.046493029	0.039476365	0.053279834	0.089581834	0.12072004	0.117061526]';
S1_krig = f_variogram_fit(input, output_S1, lb, ub);
% S2
output_S2 = [0.648459753	0.651914207	0.654000258	0.655060539	0.65588954	0.652098685	0.597778826	0.581792998	0.568663043	0.556537124	0.54366741	0.529771604	0.552524558	0.534963717	0.514637534	0.487539448	0.45385543	0.41242543	0.535606239	0.507860637	0.467424336	0.411117234	0.335185098	0.254863405	0.514618319	0.458210909	0.352781764	0.227665592	0.273764685	0.478036386]';
S2_krig = f_variogram_fit(input, output_S2, lb, ub);
% S3
output_S3 = [0.165704023	0.152346709	0.142202819	0.134298968	0.127747315	0.095685059	0.039620931	0.03399165	0.031802498	0.031148058	0.031294479	0.031928615	0.025404173	0.025119848	0.025622656	0.026656967	0.028090582	0.029864917	0.015878777	0.017832403	0.01992623	0.022675186	0.023801415	0.018764383	0.013355108	0.014182248	0.014517613	0.011746992	0.007176491	0.006431905]';
S3_krig = f_variogram_fit(input, output_S3, lb, ub);
% S4
output_S4 = [0.000571683	0.000606639	0.00064117	0.000638427	0.000563566	0.000558312	0.000463542	0.00046579	0.000433613	0.000426506	0.000414073	0.000409397	0.000426417	0.000407519	0.000384213	0.000375174	0.000390495	0.000430282	0.000364472	0.000389586	0.000426708	0.000532271	0.000634561	0.000723335	0.000406404	0.000496673	0.000657132	0.000853645	0.000623567	0.000563607]';
S4_krig = f_variogram_fit(input, output_S4, lb, ub);
% S5
output_S5 = [0.121338062	0.129753143	0.137416971	0.143588287	0.149160085	0.17969085	0.273655808	0.29857909	0.31947094	0.33967418	0.360036451	0.381007916	0.339368381	0.36805173	0.398132357	0.431396629	0.465059747	0.496519655	0.371972571	0.411098505	0.454893307	0.49785124	0.544382813	0.587070689	0.403696174	0.463831513	0.55234338	0.641873723	0.573560682	0.378583576]';
S5_krig = f_variogram_fit(input, output_S5, lb, ub);
% S6
output_S6 = [0.000238897	0.000243456	0.000295111	0.000285241	0.000236964	0.000240086	0.000225963	0.000161558	0.000204426	0.000234977	0.00023542	0.000258459	0.000199767	0.0002298	0.000232572	0.000224147	0.000255811	0.00028274	0.000209278	0.000210283	0.000244104	0.00028855	0.000284062	0.000313434	0.000200611	0.000249678	0.000294366	0.000386412	0.000484405	0.000209787]';
S6_krig = f_variogram_fit(input, output_S6, lb, ub);
% St1
output_St1 = [0.04443894	0.047321668	0.049232451	0.050461768	0.051231488	0.05943501	0.082148942	0.079717682	0.074335931	0.067077249	0.05926567	0.051686926	0.077769101	0.067012064	0.056823776	0.049570645	0.048120996	0.056152927	0.072591598	0.059268124	0.05358712	0.064063846	0.092268485	0.135150182	0.064806783	0.060145991	0.076645989	0.114985118	0.142299795	0.133947455]';
St1_krig = f_variogram_fit(input, output_St1, lb, ub);
% St2
output_St2 = [0.674287276	0.676991085	0.677938138	0.678624882	0.678877435	0.673387915	0.61555235	0.599355923	0.586649148	0.574940433	0.562746073	0.549331986	0.57071664	0.55394625	0.53450687	0.508739556	0.476183491	0.43592387	0.555173855	0.529000058	0.49016075	0.435713906	0.360669194	0.27932515	0.535780092	0.481409603	0.377891191	0.253180428	0.291803078	0.490386494]';
St2_krig = f_variogram_fit(input, output_St2, lb, ub);
% St3
output_St3 = [0.185952401	0.171257563	0.159900046	0.151183535	0.143961386	0.108613503	0.046287393	0.039898835	0.037322714	0.036502041	0.036597307	0.037253142	0.030087934	0.029631276	0.030122485	0.031212908	0.032732	0.034651863	0.019422973	0.021456969	0.023687805	0.026722502	0.027973501	0.022662971	0.016517207	0.017383059	0.017745986	0.015062184	0.010562984	0.009442522]';
St3_krig = f_variogram_fit(input, output_St3, lb, ub);
% St4
output_St4 = [0.002770566	0.002767784	0.002756615	0.002741247	0.002705826	0.002675279	0.002463113	0.002394694	0.00234635	0.002345234	0.002319337	0.00230801	0.002300183	0.002273844	0.002253021	0.002279258	0.002311363	0.002358276	0.002259736	0.002278984	0.002338044	0.002474367	0.002598942	0.002777991	0.002295881	0.002421162	0.00259987	0.002961186	0.00308874	0.002740381]';
St4_krig = f_variogram_fit(input, output_St4, lb, ub);
% St5
output_St5 = [0.123939652	0.132129885	0.139755401	0.145978609	0.151631434	0.182159041	0.275909272	0.300919649	0.321802677	0.341967925	0.362454657	0.383377256	0.341624601	0.370396599	0.400477278	0.433694504	0.467437287	0.498978213	0.374327797	0.413364381	0.457352221	0.500282501	0.54717148	0.590768645	0.405921131	0.466305239	0.555510917	0.646329691	0.580858972	0.387804206]';
St5_krig = f_variogram_fit(input, output_St5, lb, ub);
% St6
output_St6 = [0.002785233	0.002741068	0.002743121	0.002735503	0.002714117	0.002650692	0.002448373	0.002430263	0.00241442	0.002405758	0.002390966	0.002389599	0.002348276	0.002330569	0.002324836	0.002346702	0.002350069	0.002331968	0.002305862	0.002333803	0.002339052	0.002366909	0.002375278	0.002339435	0.002323546	0.002339904	0.002390817	0.002552176	0.002790439	0.00235259]';
St6_krig = f_variogram_fit(input, output_St6, lb, ub);
% Var
output_var = [37.83793989	39.85345757	16.34687101	19.00411407	21.89232913	20.01431635	10.75508513	13.38128907	16.31967576	19.41982328	22.54911806	25.62308533	15.69193347	19.96092729	24.14134391	28.06911034	31.76876424	35.47127013	21.65054202	27.14714113	31.9477636	36.27658479	41.31340392	47.37118862	26.96966126	32.3649755	36.66753059	42.93577673	68.60742236	145.6357425]';
Var_krig = f_variogram_fit(input, output_var, lb, ub);

%% Historical Data
S_1 = S1_est;
S_2 = S2_est;
S_3 = S3_est;
S_4 = S4_est;
S_5 = S5_est;
S_6 = S6_est;
St_1 = St1_est;
St_2 = St2_est;
St_3 = St3_est;
St_4 = St4_est;
St_5 = St5_est;
St_6 = St6_est;
Var = var_est;




%% plot kriging
% resolution = 100;
% [xp1, xp2] = meshgrid(lb(1,1):abs(lb(1,1)-ub(1,1))/resolution:ub(1,1), ...
%     lb(1,2):abs(lb(1,2)-ub(1,2))/resolution:ub(1,2));
% 
% xp(:,:,1) = xp1;
% xp(:,:,2) = xp2;
% xp = reshape(xp, [], 2);
% fp_krig = f_predictkrige(xp, S1_krig);
% fp_krig = reshape(fp_krig, resolution+1, resolution+1);
% 
% % figure(1)
% % hold on
% % contour(xp1, xp2, fp_krig,1000);
% % hold off
% % axis square
% figure(2)
% surf(xp1, xp2, fp_krig);
% xlabel('$x_1$','interpreter','latex');
% ylabel('$x_2$','interpreter','latex');
% 
% tic
% setting
options.conTol = 1e-6;
options.display = 0;
fileInfo.fName = 'obj';
% algorithm
result = UMDIRECT(fileInfo,lb,ub, options);
% result
x_best = result.xBest;
y_best = result.fBest;
toc

% best set
S1_best = f_predictkrige(x_best, S1_krig);
S2_best = f_predictkrige(x_best, S2_krig);
S3_best = f_predictkrige(x_best, S3_krig);
S4_best = f_predictkrige(x_best, S4_krig);
S5_best = f_predictkrige(x_best, S5_krig);
S6_best = f_predictkrige(x_best, S6_krig);
St1_best = f_predictkrige(x_best, St1_krig);
St2_best = f_predictkrige(x_best, St2_krig);
St3_best = f_predictkrige(x_best, St3_krig);
St4_best = f_predictkrige(x_best, St4_krig);
St5_best = f_predictkrige(x_best, St5_krig);
St6_best = f_predictkrige(x_best, St6_krig);
Var_best = f_predictkrige(x_best, Var_krig);

opt_S = [S1_best S2_best S3_best S4_best S5_best S6_best St1_best St2_best St3_best St4_best St5_best St6_best Var_best];

history = result.rect.x;


subplot(6,1,1)
Sobol_indice = [S1_best St1_best;S2_best St2_best;S3_best St3_best;S4_best St4_best;S5_best St5_best;S6_best St6_best];
bar(Sobol_indice);
set(gca, 'xticklabel', {'\theta_1','\theta_2','\theta_3','\theta_4','\theta_5','\theta_6','interpreter','latex','Fontsize',18});
ylabel('Sensitivity','Interpreter','latex','Fontsize',11)
xlabel('Uncertain Parameters','Interpreter','latex','Fontsize',11)
% h = '[$x_1^*,x_2^*$]=['+'num2str(x_best(1))'+', '+'num2str(x_best(2))+]';
h = ['[$x_1^*,x_2^*$]=[',num2str(x_best(1)),', ',num2str(x_best(2)),']',', Var = ',num2str(Var_best),' on exciting $\theta_1$'];
% title('[$x_1^*,x_2^*$]=[5.0001,0.1140,9.9997]','Interpreter','Latex','Fontsize',13);
title(h,'Interpreter','Latex','Fontsize',13);
legend('MSI','TSI','Interpreter','latex','Fontsize',13)
axis([0.514285714285714,6.485714285714286,0,0.65])


