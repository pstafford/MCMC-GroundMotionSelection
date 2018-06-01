% this function is used to predict mean and std of Ds 595
% reference Bommer, J.J., Stafford, P.J. and Alarcón, J.E., 2009. Empirical equations for the 
% prediction of the significant, bracketed, and uniform duration of earthquake ground motion. 
% Bulletin of the Seismological Society of America, 99(6), pp.3217-3233.
function [Ds,sigma_Ds]=BSA_09_Ds (magnitude,rupture,V30input,Ztor_input,arb_input)
R=rupture;
M=magnitude;
Ztor=Ztor_input;
Vs30=V30input;
i=2;

c0=[-5.6298,-2.2393];
m1=[1.2619,0.9368];
r1=[2.0063,1.5686];
r2=[-0.252,-0.1953];
h1=[2.3316,2.5];
v1=[-0.29,-0.3478];
z1=[-0.0522,-0.0365];

sigma_arb=[0.5564,0.4748];
sigma_geo=[0.5289,0.4616];

Ds=exp(c0(i)+m1(i)*M+(r1(i)+r2(i)*M)*log(sqrt(R^2+h1(i)^2))+v1(i)*log(Vs30)+z1(i)*Ztor);

if arb_input==1
	sigma_Ds=sigma_arb(i);
else
	sigma_Ds=sigma_geo(i);
end
