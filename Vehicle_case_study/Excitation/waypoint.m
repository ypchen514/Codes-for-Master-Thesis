function waypoint = waypoint(r)
t = linspace(0,10,600);
x = r*cos(10*t-pi/2);
y = r*sin(10*t-pi/2);
x_ini = linspace(0,r,10);
y_ini = zeros(1,length(x_ini));
x = [x_ini x+r];
y = [y_ini y+r];
waypoint = [x' y'];
plot(waypoint(:,1),waypoint(:,2),'r.');