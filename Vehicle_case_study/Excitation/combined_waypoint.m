function [waypoint, require_velocity] = combined_waypoint(sys_par)

%% Double Lane Change
sys_par = [5.8253,1.1578,0,1,1]; % sys_par for DLane
require_velocity = 1;
a1 = sys_par(1);
a2 = sys_par(2);
b = sys_par(3);
zoom_rate = sys_par(4);
zoom_long = sys_par(5);
%zoom_rate should be <=zoom_long
a1 = a1/10;
a2 = (10-a2)/10;
b = b/10;
store_data = [];
t = linspace(0,1,99);
pt1 = [0 0;a1 b;a2 (1-b);1 1]*5;
pt2 = [0 1;1-a2 1-b;1-a1 b;1 0]*5;
for ii = 1:2
    eval(['pt=pt',num2str(ii),';']);
    pts = kron((1-t).^3,pt(1,:)') + kron(3*(1-t).^2.*t,pt(2,:)') + kron(3*(1-t).*t.^2,pt(3,:)') + kron(t.^3,pt(4,:)');
    pts(1:2,end+1)=pts(1:2,end)+pts(1:2,end)-pts(1:2,end-1);
    ptsX(ii,:) = pts(1,1:length(pts));
    ptsY(ii,:) = pts(2,1:length(pts));
    store_data(2*ii-1,:) = pts(1,:);
    store_data(2*ii,:) = pts(2,:);
end
% end
points = [store_data(1,:) linspace(store_data(1,end),store_data(1,end)+store_data(3,1)+2,20)  store_data(3,:)+store_data(1,end)+2;
    store_data(2,:) ones(1,20)*store_data(2,end) store_data(4,:)];
over_dis = 3;
x_points = [linspace(0,over_dis,50) points(1,:)+over_dis linspace(points(1,end)+over_dis,points(1,end)+10*over_dis,500)]*zoom_long;
y_points = [zeros(1,50) points(2,:) ones(1,500)*points(2,end)]*zoom_rate;
waypoint_dlane = [x_points(1:300)' y_points(1:300)'];

%% Circle
sys_par = [2.0011,1]; % sys_par of circle
r = sys_par(1);
t = linspace(0,3.8,1000);
x = r*cos(2*r*t-pi/2);
y = r*sin(2*r*t-pi/2);
x_ini = linspace(0,r,10);
y_ini = zeros(1,length(x_ini));
x = [x_ini x+r];
y = [y_ini y+r];
waypoint_circle = [x'+waypoint_dlane(end,1) y'+waypoint_dlane(end,2)];
require_velocity = sys_par(2);



%% chirp
sys_par = [0.2434,0.9905,0.01];
w = sys_par(1);
A = sys_par(2);
gr = sys_par(3);
t = linspace(0,90,10000);
y = A*sin((0.5+0.05*exp(gr.*t))*w.*t);
waypoint_chirp = [-(t*0.2)'+waypoint_circle(end,1) y'+waypoint_circle(end,2)];

%% spline
t = linspace(-pi/2,pi/2,100);
x = -2*cos(t);
y = 2*sin(t)+2;

waypoint_spline = [x'+waypoint_chirp(end,1) y'+waypoint_chirp(end,2)];

%% inverse chirp
require_velocity = 1;
sim_time = 30;
% Needed Parameters: angular frequency omega,amplitude A,growing
% rate gr
sys_par = [0.2,0.7327,0.1];
w = sys_par(1);
A = sys_par(2);
gr = sys_par(3);
t = linspace(0,120,1500);
y = A*sin((1-0.00000002*exp(gr.*t))*w.*t);
t = (t*0.2)'+0.5;
tk = linspace(0,0.5,100);
t = [tk';t];
y = [zeros(1,100) y];

waypoint_invchirp = [t+waypoint_spline(end,1), y'+waypoint_spline(end,2)];



waypoint_x = [waypoint_dlane(:,1);waypoint_circle(:,1);waypoint_chirp(:,1);waypoint_spline(:,1);waypoint_invchirp(:,1)];
waypoint_y = [waypoint_dlane(:,2);waypoint_circle(:,2);waypoint_chirp(:,2);waypoint_spline(:,2);waypoint_invchirp(:,2)];

waypoint = [waypoint_x, waypoint_y];
require_velocity = 1;

plot(waypoint_x,waypoint_y)
hold on
plot(waypoint_x(end),waypoint_y(end),'r.','markersize',12)
plot(waypoint_x(1:length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)+length(waypoint_spline)),waypoint_y(1:length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)+length(waypoint_spline)));
plot(waypoint_x(1:length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)),waypoint_y(1:length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)));
plot(waypoint_x(1:length(waypoint_dlane)+length(waypoint_circle)),waypoint_y(1:length(waypoint_dlane)+length(waypoint_circle)));
plot(waypoint_x(1:length(waypoint_dlane)),waypoint_y(1:length(waypoint_dlane)));
plot(waypoint_x(length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)+length(waypoint_spline)),waypoint_y(length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)+length(waypoint_spline)),'r.','markersize',12);
plot(waypoint_x(length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)),waypoint_y(length(waypoint_dlane)+length(waypoint_circle)+length(waypoint_chirp)),'r.','markersize',12);
plot(waypoint_x(length(waypoint_dlane)+length(waypoint_circle)),waypoint_y(length(waypoint_dlane)+length(waypoint_circle)),'r.','markersize',12);
plot(waypoint_x(length(waypoint_dlane)),waypoint_y(length(waypoint_dlane)),'r.','markersize',12);