% DESCRIPTION
% Nonlinear simulation model of a notional utility rotorcraft 
% representative of UH-60 helicopter. The mathematical model is partly 
% taken from Howlett, J. J., "UH-60A Black Hawk Engineering Simulation 
% Program. Volume 1: Mathematical Model," TR 198423, NASA, Washington, DC,
% 1981. The trim, linearization routines, and overall architecture of the
% simulation follows the teachings of Dr. Joe Horn at Penn State. 
% 
% AUTHOR
% Umberto Saetti
% Assistant Professor
% Department of Aerospace Engineering 
% University of Maryland 
% saetti@umd.edu
% 
% DATE
% 4/7/2023

clear all
close all
clc

%--------------------------------------------------------------------------

% STATES             | INDICES     | DESCRIPTION
% _________________________________________________________________________
% u v w              |  1  2  3    | body velocities [ft/s]
% p q r              |  4  5  6    | angula rates [rad]
% phi theta psi      |  7  8  9    | Euler angles [rad]
% x  y  z            | 10 11 12    | position [ft]
% b0 b1c b1s b0d     | 13 14 15 16 | flapping angles [rad]
% db0 db1c db1s db0d | 17 18 19 20 | flapping angles derivatives [rad/s]
% z0 z1c z1s z0d     | 21 22 23 24 | lead-lag angles [rad]
% dz0 dz1c dz1s dz0d | 25 26 27 28 | lead-lag angles derivative [rad/s]
% lam0 lam1c lam1s   | 29 30 31    | main rotor inflow [rad]
% psi                | 32          | rotor azimuth [rad]
% dynamic twist      | 33          | dynamic twist force [lb]
% lam0T              | 34          | tail rotor inflow [rad]
% 
% CONTROL INPUTS     | INDICES     | DESCRIPTION
% _________________________________________________________________________
% delta_lat          | 1           | lateral stick 
% delta_lon          | 2           | longitudinal stick 
% delta_col          | 3           | collective stick 
% delta_ped          | 4           | pedals
% 
% OUTPUTS            | INDICES     | DESCRIPTION
% _________________________________________________________________________
% a_x a_y a_z        |  1  2  3    | body accelerations [ft/s^2]
% Q                  |  4          | main rotor torque [lb-ft^2]

%---------------------------- INITIALIZE MODEL ----------------------------

% desired flight speed in NED fame [kts -> ft/s]
vxdes=1*1.68781;
vydes=0;
vzdes=0;
% initial position in NED frame [ft]
xdes=0;
ydes=0;
zdes=0;
% desired heading [rad]
psides=0;
% desired azimuth angle of reference blade [rad] 
azdes=0;
% trim rotorcraft
[const,dt,x0,u0,Cmd0,ctrl0,t,u]=SimInit(...
    vxdes,vydes,vzdes,xdes,ydes,zdes,psides,azdes);

% linearize simulation 
[A,B]=LinSim('SimModel',x0,u0,const); 
% actuator inputs 
[controls0]=mixing(u0,const);
% trim state 
fprintf('\nABSPLUTE VELOCITY\n')
fprintf('V     = %2.2f [kts]\n',sqrt(x0(1)^2+x0(2)^2+x0(3)^2)*...
    const.fps2kts)
fprintf('\nBODY VELOCITIES\n')
fprintf('u     = %2.2f [ft/s]\n',x0(1))
fprintf('v     = %2.2f [ft/s]\n',x0(2))
fprintf('w     = %2.2f [ft/s]\n',x0(3))
fprintf('\nATTITUDE\n')
fprintf('phi   = %2.2f [deg]\n',x0(7)*const.R2D)
fprintf('theta = %2.2f [deg]\n',x0(8)*const.R2D)
fprintf('psi   = %2.2f [deg]\n',x0(9)*const.R2D)
fprintf('\nROTOR FLAPPING ANGLES\n')
fprintf('beta_0  = %2.2f [deg]\n',x0(13)*const.R2D)
fprintf('beta_1c = %2.2f [deg]\n',x0(14)*const.R2D)
fprintf('beta_1s = %2.2f [deg]\n',x0(15)*const.R2D)
fprintf('beta_0D = %2.2f [deg]\n',x0(16)*const.R2D)
fprintf('\nROTOR LEAD-LAG ANGLES\n')
fprintf('zeta_0  = %2.2f [deg]\n',x0(21)*const.R2D)
fprintf('zeta_1c = %2.2f [deg]\n',x0(22)*const.R2D)
fprintf('zeta_1s = %2.2f [deg]\n',x0(23)*const.R2D)
fprintf('zeta_0D = %2.2f [deg]\n',x0(24)*const.R2D)
fprintf('\nPILOT STICKS\n')
fprintf('LAT   = %2.2f [%%]\n',u0(1))
fprintf('LON   = %2.2f [%%]\n',u0(2))
fprintf('COL   = %2.2f [%%]\n',u0(3))
fprintf('PED   = %2.2f [%%]\n',u0(4))
fprintf('\nACTUATORS\n')
fprintf('theta_1c = %2.2f [deg]\n',controls0(1))
fprintf('theta_1s = %2.2f [deg]\n',controls0(2))
fprintf('theta_0  = %2.2f [deg]\n',controls0(3))
fprintf('theta_0T = %2.2f [deg]\n',controls0(4))

%% --------------------------- 8-STATE MODEL ------------------------------

% model-order reduction 
% retain [u v w p q r phi theta]
slow=[1:8];
fast=[13:31 34];
% reduce A matrix
As=A(slow,slow);
Asf=A(slow,fast);
Af=A(fast,fast);
Afs=A(fast,slow);
Ared8=As-Asf*inv(Af)*Afs;
% reduce B matrix
Bs=B(slow,:);
Bf=B(fast,:);
Bred8=Bs-Asf*inv(Af)*Bf;
% eigenvalues 
figure 
grid on 
hold on 
% plot(eig(A([1:8 13:20 29:31 34],[1:8 13:20 29:31 34])),'kx')
plot(eig(A),'kx')
plot(eig(Ared8),'ro')

%% --------------------------- 10-STATE MODEL -----------------------------

% model-order reduction 
% retain [u v w p q r phi theta beta_1s beta_1c]
slow=[1:8 15:16];
fast=[13:14 17:31 34];
% fast=[13:20 29:31 34];
% reduce A matrix
As=A(slow,slow);
Asf=A(slow,fast);
Af=A(fast,fast);
Afs=A(fast,slow);
Ared10=As-Asf*inv(Af)*Afs;
% reduce B matrix
Bs=B(slow,:);
Bf=B(fast,:);
Bred10=Bs-Asf*inv(Af)*Bf;
% eigenvalues 
figure 
grid on 
hold on 
plot(eig(A),'kx')
plot(eig(Ared10),'ro')
plot(eig(Ared8),'b<')
% time delay 
tauf_1c=-1/Ared10(9,9)
tauf_1s=-1/Ared10(10,10)

%% ---------------------- DECOUPLED 4-STATE MODELS ------------------------

% decoupled longitudinal dynamics 
Alon=Ared8([1 3 5 8],[1 3 5 8]);
Blon=Bred8([1 3 5 8],[2 3]);
% decoupled lateral dynamics 
Alat=Ared8([2 4 6 7],[2 4 6 7]);
Blat=Bred8([2 4 6 7],[1 4]);

%% ---------------------- DECOUPLED 5-STATE MODELS ------------------------

% decoupled longitudinal dynamics 
Alon5=Ared10([1 3 5 8 9],[1 3 5 8 9]);
Blon5=Bred10([1 3 5 8 9],[2 3]);
% decoupled lateral dynamics 
Alat5=Ared10([2 4 6 7 10],[2 4 6 7 10]);
Blat5=Bred10([2 4 6 7 10],[1 4]);

%% ------------------------ EIGENVECTOR ANALYSIS --------------------------

% scaling matrix 
S=[ones(3)          ones(3,5)*pi/180
   ones(5,3)*180/pi ones(5)];
Ared8scaled=Ared8.*S; 
% eigenvectors and eigenvalues 
[eigv,eigs]=eig(Ared8scaled); 
% eigenvectors
eigv
% eigenvlalues
diag(eigs)
% eigenvector magnitude
abs(eigv)
% eigenvector phase 
angle(eigv)*180/pi






















