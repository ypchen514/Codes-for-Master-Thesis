%% Call plot_S

tic
%% Optimized x1
load('S1_opt_data');
parfor i = 1: length(result.rect.x)
    x = result.rect.x(i,:);
    [S1(i),S2(i),S3(i),S4(i),St1(i),St2(i),St3(i),St4(i)] = plot_S(x);
end
toc 
St_Si = St1-S1;
SS = S2+S3+S4;

load('S3_opt_data');
parfor i = 1: length(result.rect.x)
    x = result.rect.x(i,:);
    [S13(i),S23(i),S33(i),S43(i),St13(i),St23(i),St33(i),St43(i)] = plot_S(x);
end
toc 
St_Si_3 = St33-S33;
SS_3 = S23+S13+S43;

subplot(4,2,1)
S1_plot_data = plot(S1,St_Si,'cyan.','Markersize',12);
axis([0,0.5,0,0.8]);
xlabel('$S_1$','Interpreter','latex');
ylabel('$S^t_1-S_1$','Interpreter','latex');
hold on
St_S1_direct_plot_data =  plot(0.3959,0.1524,'r.','Markersize',18);
St_S1_kriging_plot_data = plot(0.3606,0.1564,'ko','Markersize',7);
legend([St_S1_direct_plot_data,St_S1_kriging_plot_data],'Direct Optimization','Kriging Optimization','Interpreter','latex');
title('Interaction term of $x_2$','Interpreter','latex')
subplot(4,2,2)
S1_plot_data = plot(S1,SS,'cyan.','Markersize',12);
axis([0,0.5,0,0.8]);
xlabel('$S_1$','Interpreter','latex');
ylabel('$S_2+S_3+S_4$','Interpreter','latex');
hold on
St_S1_plot_data =  plot(0.3959,0.4235,'r.','Markersize',18);
SS1_plot_data = plot(0.3606,0.4386,'ko','Markersize',7);
legend([St_S1_plot_data,SS1_plot_data],'Direct Optimized','Kriging Optimization','Interpreter','latex');
title('Sum of Main effect exclude $x_1$','Interpreter','latex')


%% Optimized x2

subplot(4,2,3)
S1_plot_data = plot(S1,St_Si,'cyan.','Markersize',12);
axis([0,0.5,0,0.8]);
xlabel('$S_1$','Interpreter','latex');
ylabel('$S^t_1-S_1$','Interpreter','latex');
hold on
St_S1_direct_plot_data =  plot(0.3959,0.1524,'r.','Markersize',18);
St_S1_kriging_plot_data = plot(0.3606,0.1564,'ko','Markersize',7);
legend([St_S1_direct_plot_data,St_S1_kriging_plot_data],'Direct Optimization','Kriging Optimization','Interpreter','latex');
title('Interaction term of $x_1$','Interpreter','latex')
subplot(4,2,4)
S1_plot_data = plot(S1,SS,'cyan.','Markersize',12);
axis([0,0.5,0,0.8]);
xlabel('$S_2$','Interpreter','latex');
ylabel('$S_1+S_3+S_4$','Interpreter','latex');
hold on
St_S1_plot_data =  plot(0.3959,0.4235,'r.','Markersize',18);
SS1_plot_data = plot(0.3606,0.4386,'ko','Markersize',7);
legend([St_S1_plot_data,SS1_plot_data],'Direct Optimization','Kriging Optimization','Interpreter','latex');
title('Sum of Main effect exclude $x_2$','Interpreter','latex')


subplot(4,2,5)
S3_plot_data = plot(S33,St_Si_3,'cyan.','Markersize',12);
axis([0,0.6,0,0.5]);
xlabel('$S_3$','Interpreter','latex');
ylabel('$S^t_3-S_3$','Interpreter','latex');
hold on
St_S3_direct_plot_data =  plot(0.4907,0.017,'r.','Markersize',18);
St_S3_kriging_plot_data = plot(0.3941,0.0207,'ko','Markersize',7);
legend([St_S3_direct_plot_data,St_S3_kriging_plot_data],'Direct Optimization','Kriging Optimization','Interpreter','latex');
title('Interaction term of $x_3$','Interpreter','latex')

subplot(4,2,6)
S3_plot_data = plot(S33,SS_3,'cyan.','Markersize',12);
axis([0,0.6,0,0.8]);
xlabel('$S_3$','Interpreter','latex');
ylabel('$S_1+S_2+S_4$','Interpreter','latex');
hold on
St_S3_plot_data =  plot(0.4907,0.4912,'r.','Markersize',18);
SS3_plot_data = plot(0.3941,0.6349,'ko','Markersize',7);
legend([St_S3_plot_data,SS3_plot_data],'Direct Optimized','Kriging Optimization','Interpreter','latex');
title('Sum of Main effect exclude $x_3$','Interpreter','latex')


subplot(4,2,7)
S3_plot_data = plot(S33,St_Si_3,'cyan.','Markersize',12);
axis([0,0.6,0,0.5]);
xlabel('$S_4$','Interpreter','latex');
ylabel('$S^t_4-S_4$','Interpreter','latex');
hold on
St_S3_direct_plot_data =  plot(0.4907,0.017,'r.','Markersize',18);
St_S3_kriging_plot_data = plot(0.3941,0.0207,'ko','Markersize',7);
legend([St_S3_direct_plot_data,St_S3_kriging_plot_data],'Direct Optimization','Kriging Optimization','Interpreter','latex');
title('Interaction term of $x_4$','Interpreter','latex')

subplot(4,2,8)
S3_plot_data = plot(S33,SS_3,'cyan.','Markersize',12);
axis([0,0.6,0,0.8]);
xlabel('$S_4$','Interpreter','latex');
ylabel('$S_1+S_2+S_3$','Interpreter','latex');
hold on
St_S3_plot_data =  plot(0.4907,0.4912,'r.','Markersize',18);
SS3_plot_data = plot(0.3941,0.6349,'ko','Markersize',7);
legend([St_S3_plot_data,SS3_plot_data],'Direct Optimized','Kriging Optimization','Interpreter','latex');
title('Sum of Main effect exclude $x_4$','Interpreter','latex')