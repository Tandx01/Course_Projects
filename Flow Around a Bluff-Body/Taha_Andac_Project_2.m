clc, clear all, close all
% computational domain creation
xl=4;
yl =1;
nx = 200;
ny = 50;
uin = 1;
pexit = 0;
re = 20;
dt=0.01;
nstep=500;
maxit=1000;
beta=1.2;
h= xl/nx;
Perrtool =0.001;
fno = 1; %figure number
%%

u=zeros(nx+1,ny+2);v=zeros(nx+2,ny+1);p=zeros(nx+2,ny+2);
ut=zeros(nx+1,ny+2);vt=zeros(nx+2,ny+1);c=zeros(nx+2,ny+2)+0.25;
uu=zeros(nx+1,ny+1);vv=zeros(nx+1,ny+1);w=zeros(nx+1,ny+1);
T=zeros(nx+2,ny+2);T((nx/4)+1:(nx/2)+1,1:(ny/2)+1)=1;

c(2,3:ny)=1/3;
c(3:nx,2)=1/3;
c(3:nx,ny+1)=1/3;
c(2,2)=1/2;
c(2,ny+1)=1/2;
c(nx+1,2)=0.5;

un=0;us=0;ve=0;vw=0;
mu = uin*2*yl/re;
time=0.0;

F=[0];
Cd=[0];
Nu=[0];
KK=0;

%block edges
ibw =nx/4; %west edge
ibe = nx/2; %east edge
jbs = 1; % south edge
jbn = ny/2; %north edge

%blok pressure BC

c(ibw+1,2:jbn+1) = 1/3;
c(ibe+1,2:jbn+1) =1/3;
c(ibw+2:ibe+1,jbn+2) =1/3;
c(ibw+2:ibe,2:jbn+1) =0;

for is=1:nstep

u(1:nx+1,1) = u(1:nx+1,2); v(1:nx+1,1) =0; % symmetry south bc  
u(1:nx+1,ny+2) = -u(1:nx+1,ny+1); v(1:nx+1,ny+1) = 0; %no-slip north bc
u(1,1:ny+2) = uin; v(1,1:ny+1) =-v(2,1:ny+1); %inlet bc 
u(nx+1,1:ny+2) = u(nx,1:ny+2); v(nx+2,1:ny+1)=-v(nx+1,1:ny+1); %outflow

%blok

u(ibw+1,1:jbn+1)=0; v(ibw+2,1:jbn+1)=-v(ibw+1,1:jbn+1); %left boundary of the block
u(ibe+1,1:jbn+1) =0; v(ibe+1,1:jbn+1)=-v(ibe+2,1:jbn+1); %right boundary of blok
u(ibw+1:ibe+1,jbn+1) = -u(ibw+1:ibe+1,jbn+2); v(ibw+1:ibe+1,jbn+1)=0; %top of blok
p(ibw+2:ibe+1,2:jbn+1) =0;

Maxvel = max(max(uu.^2 + vv.^2));
Maxvel = max(1.e-6,Maxvel);
dtvis = 0.25*h*h/mu; dtadv = 4*mu/Maxvel; dt = min(dtvis,dtadv);

for i=2:nx, for j=2:ny+1 % temporary u-velocity
ut(i,j)=u(i,j)+dt*(-(0.25/h)*((u(i+1,j)+u(i,j))^2-(u(i,j)+...
u(i-1,j))^2+(u(i,j+1)+u(i,j))*(v(i+1,j)+...
v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+...
(mu/h^2)*(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)));
end, end
for i=2:nx+1, for j=2:ny % temporary v-velocity
vt(i,j)=v(i,j)+dt*(-(0.25/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+...
v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+...
(v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)+...
(mu/h^2)*(v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)));
end, end
% apply bc for vt and ut
% velocity BC 
ut(1:nx+1,1) = ut(1:nx+1,2); vt(1:nx+1,1) =0; % symmetry south bc  
ut(1:nx+1,ny+2) = -ut(1:nx+1,ny+1); vt(1:nx+1,ny+1) = 0; %no-slip north bc
ut(1,1:ny+2) = uin; vt(1,1:ny+1) =-vt(2,1:ny+1); %inlet bc 
ut(nx+1,1:ny+2) = ut(nx,1:ny+2); vt(nx+2,1:ny+1)=-vt(nx+1,1:ny+1); %outflow

%blok
ut(ibw+1,1:jbn+1)=0; vt(ibw+2,1:jbn+1)=- vt(ibw+1,1:jbn+1); %left boundary of the block
ut(ibe+1,1:jbn+1) =0; vt(ibe+1,1:jbn+1)=-vt(ibe+2,1:jbn+1); %right boundary of blok
ut(ibw+1:ibe+1,jbn+1) = -ut(ibw+1:ibe+1,jbn+2); vt(ibw+1:ibe+1,jbn+1)=0; %top of blok
ut(ibw+2:ibe+1,2:jbn+1) =0;
vt(ibw+2:ibe+1,2:jbn+1) =0;

%pressure solver
it =0; Perr = 1.e5;
while it <maxit && Perr>Perrtool
    p(end-1,:)=pexit;
    p(ibw+2:ibe+1,2:jbn+1) =0;
    it = it+1; pold=p;
for i=2:nx
    for j=2:ny+1
        
p(i,j)=beta*c(i,j)*(p(i+1,j)+p(i-1,j)+p(i,j+1)+p(i,j-1)-...
(h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+(1-beta)*p(i,j);
    end
end

Perr =norm(p-pold);
end
% correct the velocity
u(2:nx,2:ny+1)=...
ut(2:nx,2:ny+1)-(dt/h)*(p(3:nx+1,2:ny+1)-p(2:nx,2:ny+1));
v(2:nx+1,2:ny)=...
vt(2:nx+1,2:ny)-(dt/h)*(p(2:nx+1,3:ny+1)-p(2:nx+1,2:ny));
time=time+dt;

% Force computation
    F_West=(h*(sum(p((nx/4),1:(ny/2)+1))-0.5*p((nx/4),1)-0.5*p((nx/4),(ny/2)+1)))*1; % Forces on west of object                                   
    F_East=(h*(sum(p((nx/2)+2,1:(ny/2)+1))-0.5*p((nx/2)+2,1)-0.5*p((nx/2)+2,(ny/2)+1)))*1; % Forces on east of object 
    F_North=(mu*1*1*sum((uu((nx/4)+1:(nx/2)+1,(ny/2)+2)-uu((nx/4)+1:(nx/2)+1,(ny/2)+1))/h)); % Forces on north of object 
    F(is+1)=2*(F_West-F_East+F_North); % symterty multiplication
    Cd(is+1)=F(is+1)/(0.5*1*1*1^2);
    td(is+1) = time;
    % Energy equation
    for it=1:maxit
        for i=2:nx+1
            for j=2:ny+1
                T(i,j)=T(i,j)+dt*((mu/h^2)*(T(i+1,j)+T(i-1,j)+T(i,j+1)+T(i,j-1)-4*T(i,j))-(1/(2*h))*(uu(i,j)*(T(i+1,j)-T(i-1,j))+vv(i,j)*(T(i,j+1)-T(i,j-1))));
            end
        end
        T(1:nx+1,ny+1)=-T(1:nx+1,ny);
        T(1:nx+1,2)=-T(1:nx+1,3);
        T(nx+1,1:ny+1)=-T(nx,1:ny+1);
        T((nx/4)+1:(nx/2)+1,1:(ny/2)+1)=1;
    end
    
    Nu_North=1*(sum(T((nx/4)+1:(nx/2)+1,(ny/2)+3)-T((nx/4)+1:(nx/2)+1,(ny/2)+1)))/(2*h)/(sum(T((nx/4)+1:(nx/2)+1,(ny/2)+2)-T((nx/4)+1:(nx/2)+1,(ny/2)+1)));
    Nu_West=1*(2*sum(T((nx/4)-1,1:(ny/2)+1)-T((nx/4)+1,1:(ny/2)+1)))/(2*h)/(2*sum(T((nx/4),1:(ny/2)+1)-T((nx/4)+1,1:(ny/2)+1)));
    Nu_East=1*(2*sum(T((nx/2)+3,1:(ny/2)+1)-T((nx/2)+1,1:(ny/2)+1)))/(2*h)/(2*sum(T((nx/2)+2,1:(ny/2)+1)-T((nx/2)+1,1:(ny/2)+1)));
    Nu_South=Nu_North;
    Nu(is+1)=Nu_North+Nu_South+Nu_West+Nu_East;
 
end
%%
% plot results
uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
figure(fno); fno = fno+1;
hold on
quiver(uu',vv','r');
contour(p',60),patch([(nx/4)+1 (nx/4)+1 (nx/2)+1 (nx/2)+1],[1 (ny/2)+1 (ny/2)+1 1],"k"),axis equal
colorbar
title('Pressure contours and Velocity vectors')
%%
figure(fno); fno = fno+1;
contour(T',60),patch([(nx/4)+1 (nx/4)+1 (nx/2)+1 (nx/2)+1],[1 (ny/2)+1 (ny/2)+1 1],"k"),axis equal
colorbar
title('Temperature')
%%
figure(fno); fno = fno+1;

%Computing vorticity
w(1:nx+1,1:ny+1)=(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1)-...
v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1))/(2*h);
contour(w',60),patch([(nx/4)+1 (nx/4)+1 (nx/2)+1 (nx/2)+1],[1 (ny/2)+1 (ny/2)+1 1],"k"),axis equal
colorbar
title('Vorticiy')
%%
% Compute stream function
psi = zeros(nx+2, ny+2);
for i = 2:nx+1
    for j = 2:ny+1
        psi(i, j) = psi(i-1, j) + (u(i, j) + u(i-1, j)) * 0.5 * h;
    end
end

% Plot stream function contours
figure(fno); fno = fno + 1;
contour(psi', 60), patch([(nx/4)+1 (nx/4)+1 (nx/2)+1 (nx/2)+1], [1 (ny/2)+1 (ny/2)+1 1], "k"), axis equal
colorbar
title('Stream Function')
%%
%Force versus Time Plot
figure(fno); fno = fno + 1;
plot(td,F,'k-',LineWidth=2)
xlabel('Time (s)')
ylabel('Force (N)')
grid on