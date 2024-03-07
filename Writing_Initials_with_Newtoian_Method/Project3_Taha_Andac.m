clear all, close all, clc
disp('Writing T_A:')

tol = 0.00001;
max_it = 50;
d1 = 1.5; % lenght of blue arm
d2 = 1.5; % lenght of red arm
n = 100; %step number of sub parts
t = linspace(0,1,n); 
N = 8*n; %step total step number
T = linspace(0,21,N); % Is chosen as 21 because animation last 21 seconds.
px1 = 0.5+t; 
py1 = ones(1,n); 
px2 = 1.5-0.5*t;
py2 = ones(1,n);
px3 = ones(1,n);
py3 = 1-t;
px4 = 1+0.5*t;
py4 = zeros(1,n);
px5 = 1.5+ 0.5*t;
py5 = 0+t;
px6 = 2+0.5*t;
py6 = 1-t;
px7 = 2.5-0.25*t;
py7 = 0.5*t;
px8 = 2.25-0.5*t;
py8 = ones(1,n)*0.5;

x = [px1 px2 px3 px4 px5 px6 px7 px8]; %  x coordinate path
y = [py1 py2 py3 py4 py5 py6 py7 py8]; %  y coordinate path


alpha = zeros(1,N);
beta = zeros(1,N);
angles = [0.5 0.5]'; % Initial guess
for i = 1:1:N
iter = 1;
while (iter <= max_it)
    a = angles(1); % a = Angle alpha
    b = angles(2); % b = Angle beta
    f = [ (d1*cos(a)+d2*cos(a+b)-x(i))
    (d1*sin(a)+d2*sin(a+b)-y(i)) ];
    J = [ (-d1*sin(a)-d2*sin(a+b)) (-d2*sin(a+b))
    ( d1*cos(a)+d2*cos(a+b)) ( d2*cos(a+b)) ];
    angles_new = angles -J\f;
    err = norm(angles_new - angles);
    
    if err <= tol
        angles_next = angles_new;
        break
    else
        angles = angles_new;
    end
    iter = iter + 1;
end
if (iter > max_it), error('Newton method did not converge'); end
angles = mod(angles_next,2*pi);
alpha(i) = angles(1);
beta(i) = angles(2);
end
RESULT_t_alpha_beta = [T' alpha' beta']
x0 = zeros(n);
y0 = zeros(n);
x1 = d1*cos(alpha);
y1 = d1*sin(alpha);
x2 = x1+d2*cos(alpha+beta);
y2 = y1+d2*sin(alpha+beta);

figure(1)
i = 1;
plot([x0(i) x1(i)],[y0(i) y1(i)],'b-', [x1(i) x2(i)],[y1(i) y2(i)],'r-',x2(i),y2(i),'ko','linewidth',3), 
grid on, axis([-1.5 3.5 -0.5 2.5]), xlabel('x'), ylabel('y'), 
text_plot = ['time = ',num2str(T(1))];
text(0,-0.25,text_plot);
title('Positioning of Robot Arm at T = 0'),

% Important note, sometimes MATLAB plots one point of figure 2 on figure 1
% If that happens, I recommend re-running the code, or commenting whole
% figure 1 part.

figure(2)
for i=1:N
    
hold on
text_plot = ['time = ',num2str(T(i))];
A = plot([x0(i) x1(i)],[y0(i) y1(i)],'b-', [x1(i) x2(i)],[y1(i) y2(i)],'r-','linewidth',3);
grid on, axis([-1.5 3.5 -0.5 2.5])
B = text(0,-0.25,text_plot);

K = plot(x2(i),y2(i),'ks',LineWidth=4);
grid on, axis([-1.5 3.5 -0.5 2.5]), xlabel('x'), ylabel('y'),
title('Writing T-A'), axis([-1 3 -0.5 2.5])
if i ~= N
pause(0.25)
delete (A), delete (B)
end
end
hold off

figure(4)
i = 80;
plot([x0(i) x1(i)],[y0(i) y1(i)],'b-', [x1(i) x2(i)],[y1(i) y2(i)],'r-',x2(i),y2(i),'ko','linewidth',3), 
grid on, axis([-1.5 3.5 -0.5 2.5]), xlabel('x'), ylabel('y'),
text_plot = ['time = ',num2str(T(end))];
text(0,-0.25,text_plot);
title('Positioning of Robot Arm at T = 21'),

figure(5)
alpha_deg = alpha/(2*pi)*360; % convert alpha into degrees
beta_deg = (beta/(2*pi)*360); %covert beta into degrees
subplot(2,1,1)
plot(T,alpha_deg,'b-',LineWidth=3),grid on
xlim([0 21]), title('Alfa versus time')
subplot(2,1,2)
plot(T,beta_deg,'r-',LineWidth=3),grid on
xlim([0 21]), title('Beta versus time')
