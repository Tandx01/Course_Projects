%You can just run the code without openning the simulink file.
clc,clear all,close all

s = tf('s');

%System values
K = 1;
u = 0.05;
umin = 0.2;
J = 100;
m = 0.01;
g = 9.81;

%Retrieving KP-KD values ​​found for each system
import_file("PD")
PD_tuned = PD;
Kd = PD.K;
Kp = -cell2mat(PD.Z)*Kd;

import_file("marginaly_stable")
PD_M = marginaly_stable;
Kd_M = PD_M.K;
Kp_M = -cell2mat(PD_M.Z)*Kd;

import_file("underdamped")
PD_U = underdamped;
Kd_U = PD_U.K;
Kp_U = -cell2mat(PD_U.Z)*Kd;

import_file("overdamped")
PD_O = overdamped;
Kd_O = PD_O.K;
Kp_O = -cell2mat(PD_O.Z)*Kd;

%Running simulink file
out = sim("Part_2_Simulink.slx");

%Retrieving time and position values from the simulink model
teta = out.teta;
position = out.position;
t =linspace(0,12,length(position));

teta_M = out.teta_M;
position_M = out.position_M;

teta_U = out.teta_U;
position_U = out.position_U;

teta_O = out.teta_O;
position_O = out.position_O;

%Creation of Transfer function step-by-step
fprintf("Transfer function of DC motor, Q(s)/I(s) = \n \n")
GDc = K/(J*s^2)

fprintf("Transfer function of mechanical part, X(s)/Q(s) = \n \n")
G = m*g/(m*s^2+u*s+umin)

fprintf("Open loop transfer function without controller, X(s)/I(s) = \n \n")
L = PD_tuned*GDc*G

fprintf("Closed loop transfer function, X(s)/Xdes(s) = \n \n")
T = minreal(L/(1+L))

figure(1)
subplot(2,1,1)
step(T)
title("Step response for Kp = " +Kp+ " and Kd = " +Kd)
stepinfo(T)
subplot(2,1,2)
rlocus(T)
title("Root locus for Kp = " +Kp+ " and Kd = " +Kd)
%%
%Simulation of the ball with given function
figure(2)
drawBallnBeam(position,teta)
%%
%Marginally stable system

L = PD_M*GDc*G;
T = minreal(L/(1+L));

figure(3)
pzmap(T)
title("Marginally stable root locus for Kp = " +Kp_M+ " and Kd = " +Kd_M)

figure(4)
step(T,12)
title("Marginally stable step response for Kp = " +Kp_M+ " and Kd = " +Kd_M)

figure(5)
subplot(2,1,1)
plot(t,teta_M)
xlabel("time")
ylabel("teta")
title("Marginally stable system Kp = " +Kp_M+ " and Kd = " +Kd_M)
grid on

subplot(2,1,2)
plot(t,position_M)
xlabel("time")
ylabel("position")
grid on
%%
%Underdamped
L = PD_U*GDc*G;
T = minreal(L/(1+L));
figure(6)
pzmap(T)
title("Underdamped root locus for Kp = " +Kp_U+ " and Kd = " +Kd_U)

figure(7)
step(T,12)
title("Underdamped step response for Kp = " +Kp_U+ " and Kd = " +Kd_U)

figure(8)
subplot(2,1,1)

plot(t,teta_U)
xlabel("time")
ylabel("teta")
title("Underdamped system Kp = " +Kp_U+ " and Kd = " +Kd_U)
grid on

subplot(2,1,2)
plot(t,position_U)
xlabel("time")
ylabel("position")
grid on

%%
%Overdamped
L = PD_O*GDc*G;
T = minreal(L/(1+L));
figure(9)
pzmap(T)
title("Overdamped root locus for Kp = " +Kp_O+ " and Kd = " +Kd_O)

figure(10)
step(T,12)
title("Overdamped step response for Kp = " +Kp_O+ " and Kd = " +Kd_O)

figure(11)
subplot(2,1,1)

plot(t,teta_O)
xlabel("time")
ylabel("teta")
title("Overdamped system Kp = " +Kp_O+ " and Kd = " +Kd_O)
grid on

subplot(2,1,2)

plot(t,position_O)
xlabel("time")
ylabel("position")
grid on