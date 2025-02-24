clear all, close all, clc

dt=0.001;
PopSize=25;
MaxGenerations=10;
s = tf('s');
G= 0.397/(s*(s^2+0.9471*s+0.3943))
options = optimoptions(@ga, 'PopulationSize', PopSize, 'MaxGenerations', MaxGenerations);
[x,fval] = ga(@(K)pidtest(G,dt,K),3,-eye(3),zeros(3,1))