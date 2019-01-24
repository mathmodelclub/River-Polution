%Polluted River Project
%Date: 21.11.2018
%Description: 1D Diffusion - Transport 
%Author: Alber Lukas/ Alber David/ Tollet Markus
close all
clear all
%% Load in initial parameters,  ===========================================
L       = 500; %[km]     length of river
a       = 400; %[km]     position of leakage
da      = 75;  %[km]     width of leakage
D       = 25;  %[km^2/h] difusitivity pesticide to water
v       = 5;   %[km/h]   velocity of river
psi0    = 50;  %[mol]    amplitude of leakage
alpha   = 10;  %[1/km]   diffusitivity from river to ocean 
tswitch = 10;  %[h]      time when the source is turned off
p0      = 0.6; %[mol/km] U_infty threshold
p1      = 5;   %[mol*h/km] U_tot threshold
rho     = 1;   %[1/km]  amount of fish per length
%% Discretize space and time===============================================
x0 = 0; xend = L;   hx  = 5; 
t0 = 0; tend = 150; ht  = 0.1;
vstart = 0; vend = 10; hv = 0.1;
%create space and time grid
xx = plib.xgrid(x0,xend,hx); 
tt = plib.tgrid(t0,tend,ht); 
xsize = length(xx);
tsize = length(tt);
%set initial conditions for the problem
u0 = zeros(xsize,1);
%%
% Main evaluation =========================================================
%set the right coefficients which are pluged in Matrix
[c_m1,c_0,c_p1,omega] = plib.coeffconstr(hx,v,D,alpha);
%cunstruct the sparse matrix A for ODE
T = plib.matrixconstr(xsize,c_m1,c_0,c_p1,omega);
%construct source => psix for t<tswitch; zeros t>tswitch
[psix,source] = plib.source(xx,tt,psi0,tswitch,a,da);
%%
%Creating Matrix for CN:
[A_inv,B] = plib.CNmatrix(T,xsize,ht);
%solve with CN
[U] = plib.CN(A_inv,B,xsize,tsize,source,tt,ht,u0);
%%
%U SURFPLOT
plib.SurfPlot(xx,tt,U,'txt')
%U ANIMATION
plib.animation(xx,tt,U,0)

%Source Plot
%{
[xm,tm]=meshgrid(xx,tt);
figure
surf(xm',tm',source,'edgecolor','none')
xlabel('Space x')
ylabel('Time t')
zlabel('Concentration')
title('Source(x,t)')
%}
%ONE VELOCITY U_INFTY/ U_TOT PLOT
%{
[deadinfty,U_infty] = plib.Uinfty(U,p0,rho,xx,hx);
[deadtot,U_tot] = plib.Utot(tt,U,p1,rho,xx,hx);
plib.U_inftyPlot(xx,U_infty,p0,rho,deadinfty)
plib.U_totPlot(xx,U_tot,p1,rho,deadtot)
%}

%%
%MULTIVELOCITY: Dead Fish for different velocities
[UV] = plib.VSol(vstart,vend,hv,tt,tsize,ht,xx,xsize,hx,D,alpha,source,u0);
[deadinftyV,U_inftyV] = plib.UinftyV(vstart,vend,hv,xx,hx,UV,p0,rho);
[deadtotV,U_totV] = plib.UtotV(vstart,vend,hv,xx,hx,tt,UV,p1,rho);

%%
%Animation: Dead Fish for Total Pesticide Exposure measurement
plotdim = [0 500 0 65];
plib.UVanimation(U_totV,deadtotV,vstart,vend,hv,xx,p1,rho,plotdim)
%%
%Animation: Dead Fish for Maximal Pesticide Exposure measurement
plotdim = [0 500 0 1.5];
plib.UVanimation(U_inftyV,deadinftyV,vstart,vend,hv,xx,p0,rho,plotdim)

%%
%TIEM 2 OCEAN
[Otime] = plib.Time2Ocean(vstart,vend,hv,UV,tt);
plib.OceanvecPlt(vstart,vend,hv,Otime)










