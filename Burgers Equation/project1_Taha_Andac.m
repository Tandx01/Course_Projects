clear all, clc, close all

dt = 0.001;
t = 0:dt:1;
m = length(t); % time grid points

n = 201; % space grid points
%which means we have 200 space grid
x = linspace(0,1,n);
dx = 1/(n-1); %Space grid lenght

nu = 0.01; % viscosity

u(1,:) = 1 + sin(2*pi*x); % conservative form
unc = u; %non conservative form

%finite difference factors
fac1 = 0.5*dt/dx; 
fac2 = nu*dt/(dx*dx);
fac3 = 0.5*fac1;

    for j = 1:m-1 % time increment

    for i = 2:n-1 % space increment

        % FTCS for convervative form
        u(j+1,i) = u(j,i)-fac3*(u(j,i+1)^2 - u(j,i-1)^2) + fac2*(u(j,i+1)-2*u(j,i)+u(j,i-1));

        % FTCS for non conservative form but it diverges so was not used
        unc(j+1,i) = unc(j,i)-fac1*unc(j,i)*(unc(j,i+1)-unc(j,i-1))+fac2*(unc(j,i+1)-2*unc(j,i)+unc(j,i-1));

    end
        %Applying periodic boundary conditions
        u(j+1,1) = u(j+1,n-1);
        u(j+1,n) = u(j+1,2);
        unc(j+1,1) = unc(j+1,n-1);
        unc(j+1,n) = unc(j+1,2);
        %ploting

        if mod (j, 100) == 0
        h = plot(x,u(j,:),'k-',x,unc(j,:),'m--',linewidth = 2);
        xlabel('Spatial co-ordinate (x) \rightarrow')
        ylabel('Transport property profile (u) \rightarrow')
            title({['1-D Burgers'' equation (\nu = ',num2str(nu),')'];['time(\itt) = ',num2str(dt*j)];['Space grid = ',num2str(n-1)]})
        axis([0 1 0 2])
        grid on
        legend('Conservative Solution','Non-Conservative Solution')
        drawnow;
        refreshdata(h)
        end
    end
%%
%Calculating integrated difference for t = 1
clear all, clc, close all

x0 = 0.5; % error calculation point 

dt = 0.001;
tfinal = 1;
t = 0:dt:1;
m = length(t); % time grid points

nu = 0.01; % viscosity

res = [20:10:200]+1;

for k = 1:length(res)
    
    clear u unc

    n = res(k);

    x = linspace(0,1,n);
    dx = 1/(n-1); %Space grid lenght

    %finite difference factors
    fac1 = 0.5*dt/dx;
    fac2 = nu*dt/(dx*dx);
    fac3 = 0.5*fac1;

    u(1,:) = 1 + sin(2*pi*x); % conservative form
    unc = u; %non conservative form

    i0 = find(x0 ==x); % error point index

    for j = 1:m-1 % time increment

    for i = 2:n-1 % space increment

        % FTCS for convervative form
        u(j+1,i) = u(j,i)-fac3*(u(j,i+1)^2 - u(j,i-1)^2) + fac2*(u(j,i+1)-2*u(j,i)+u(j,i-1));

        % FTCS for non conservative form but it diverges so was not used
        unc(j+1,i) = unc(j,i)-fac1*unc(j,i)*(unc(j,i+1)-unc(j,i-1))+fac2*(unc(j,i+1)-2*unc(j,i)+unc(j,i-1));

    end
        %Applying periodic boundary conditions
        u(j+1,1) = u(j+1,n-1);
        u(j+1,n) = u(j+1,2);
        unc(j+1,1) = unc(j+1,n-1);
        unc(j+1,n) = unc(j+1,2);

    end

a(k) = u(end, i0); %dumy variable for storing conservative solution
anc(k) = unc(end, i0); %dumy variable for storing non-conservative solution

end
res = 1./res;

figure(2)
plot(res,a,'k*-', res,anc,'mo--',LineWidth=2)
xlabel ('h')
ylabel('U value')
title('Accuracy at point x =',num2str(x0))
legend('Conservative Solution at different grids','Non-Conservative Solution at different grids')
grid on