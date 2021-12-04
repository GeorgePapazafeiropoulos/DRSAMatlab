%% 3D one-story building subjected to the El Centro earthquake motion (Chopra, 2012)

%% Statement of the problem
% * *Chopra (2012), Example 10.6:* Figure 9.5.1 shows the plan view of a
% one-story building. The structure consists of a roof, idealized as a
% rigid diaphragm, supported on three frames, A, B, and C, as shown. The
% roof weight is uniformly distributed and has a magnitude of $$100
% lb/ft.^2$. The lateral stiffnesses of the frames are $$k_y = 75 kips/ft$
% for frame A, and $$k_x = 40 kips/ft$ for frames B and C. The plan
% dimensions are $$b = 30 ft$ and $$d = 20 ft$, the eccentricity is $$e =
% 1.5 ft$, and the height of the building is $$12 ft$. Determine the
% natural periods and modes of vibration of the structure.
% * *Chopra (2012), Example 13.12:* Determine the peak values of the
% response of the one-story unsymmetric-plan system of Examples 13.7 and
% 10.6 with modal damping ratios $$\mathrm{\zeta_n} = 0.05$ to the El
% Centro ground motion in the y-direction, directly from the response
% spectrum for this ground motion.
%
% <<Chopra951.png>>
%
%% Initialization of structural input data
% Set the dimensions of the rigid slab (ft).
b=30;
d=20;
%%
% Set the height of the building (ft).
h=12;
%%
% Set the eccentricity (ft).
e=1.5;
%%
% Set the lateral stiffnesses of the frames (kips/ft).
ky=75;
kx=40;
%%
% Set the mass of the slab (kips*sec^2/ft).
m=1.863;
%%
% Set the moment of inertia of the roof about its center of mass
% (kips*ft*sec^2).
Io=m*(b^2+d^2)/12;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure according to eq.(9.5.17).
K=[2*kx,0,0;0,ky,e*ky;0,e*ky,e^2*ky+(d^2/2)*kx];
%%
% Calculate the mass matrix of the structure.
M=[m,0,0;0,m,0;0,0,Io];
%%
% Calculate the length of the semi-diagonal.
rbd=sqrt(b^2+d^2)/2;
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied only at dof No 2 of the structure.
r=[0;1;0];
%% Load earthquake response spectrum data
% Earthquake displacement response spectrum of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
Sd=load('Sd3.txt');
%%
% Earthquake pseudoacceleration response spectrum of the El Centro
% earthquake (El Centro, 1940, El Centro Terminal Substation Building)
PSa=load('PSa3.txt');
%%
% Set the critical damping ratio of the response spectra
% ($$\mathrm{\xi}=0.05$)
ksi=0.05;
%% Dynamic Response Spectrum Analysis (DRSA)
% Consider all eigenmodes of the building
eigInd=[1;2;3];
%%
% Perform DRSA analysis
[U,~,~,f,omega,Eigvec] = DRSA(K,M,r,PSa,Sd,[],[],eigInd);
%%
% Set the number of eigenmodes of the structure.
neig=numel(omega);
%%
% Plot the natural modes of vibration of the structure
FigHandle=figure('Name','Natural Modes','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;b;b;0;0],[0;0;d;d;0],'LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 0 1],'markeredgecolor','k');
    hold on
    phi=10*Eigvec(3,i);
    plot([10*Eigvec(1,i)+(rbd*sin(phi)*cos(pi/4+phi/2)/cos(phi/2));
        b+10*Eigvec(1,i)+(rbd*sin(phi)*cos(pi/4+phi/2)/cos(phi/2));
        b+10*Eigvec(1,i)-(rbd*sin(phi)*cos(pi/4+phi/2)/cos(phi/2));
        10*Eigvec(1,i)-(rbd*sin(phi)*cos(pi/4+phi/2)/cos(phi/2));
        10*Eigvec(1,i)+(rbd*sin(phi)*cos(pi/4+phi/2)/cos(phi/2))],...
        [10*Eigvec(2,i)-(rbd*sin(phi)*sin(pi/4+phi/2)/cos(phi/2));
        10*Eigvec(2,i)+(rbd*sin(phi)*sin(pi/4+phi/2)/cos(phi/2));
        d+10*Eigvec(2,i)+(rbd*sin(phi)*sin(pi/4+phi/2)/cos(phi/2));
        d+10*Eigvec(2,i)-(rbd*sin(phi)*sin(pi/4+phi/2)/cos(phi/2));
        10*Eigvec(2,i)-(rbd*sin(phi)*sin(pi/4+phi/2)/cos(phi/2))],...
        'LineWidth',2.,'Marker','.','MarkerSize',20,'Color',[1 0 0],...
        'markeredgecolor','k');
    hold off
    axis equal
    grid on
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Compare the eigenmodes with those shown in Figure E10.6. in Chopra (2012).
%
% <<ChopraE106.png>>
%
%%
% Convert peak modal displacements from ft to inch. Only displacements
% are multiplied. The rotation angles are given in rad.
U(1:2,:)=12*U(1:2,:);
%%
% Calculate the peak modal base shear of frame A (kips).
VbA=zeros(neig,1);
for i=1:neig
    VbA(i)=[0 ky e*ky]*U(:,i);
end
VbA
%%
% Calculate the peak modal base shear of frame B (kips).
VbB=zeros(neig,1);
for i=1:neig
    VbB(i)=[kx 0 -d/2*kx]*U(:,i);
end
VbB
%% Verification of table E13.12a of Chopra (Dynamics of Structures, 2012)
% Compare the peak modal responses with those shown in table E13.12a in
% Chopra (2012). Mode 1 corresponds to translation in y direction (here as
% mode 1) and mode 2 corresponds to rotation about z axis (here as mode 3).
%
% <<ChopraE1312a.png>>
%
%% Modal combination with the ABSolute SUM (ABSSUM) method.
% Peak displacement in y direction (inch).
uyAbsSum=ABSSUM(U(2,:)');
%%
% Peak torsional rotation (rad).
uthAbsSum=ABSSUM(U(3,:)');
%%
% Peak base shear of the whole structure (kips).
VbAbsSum=ABSSUM(f(2,:)');
%%
% Peak base torsional moment of the whole structure (kips*ft).
TbAbsSum=ABSSUM(f(3,:)');
%%
% Peak base shear of frame A (kips).
VbAAbsSum=ABSSUM(VbA);
%%
% Peak base shear of frame B (kips).
VbBAbsSum=ABSSUM(VbB);
%% Modal combination with the Square Root of Sum of Squares (SRSS) method.
% Peak displacement in y direction (inch).
uySRSS=SRSS(U(2,:)');
%%
% Peak torsional rotation (rad).
uthSRSS=SRSS(U(3,:)');
%%
% Peak base shear of the whole structure (kips).
VbSRSS=SRSS(f(2,:)');
%%
% Peak base torsional moment of the whole structure (kips*ft).
TbSRSS=SRSS(f(3,:)');
%%
% Peak base shear of frame A (kips).
VbASRSS=SRSS(VbA);
%%
% Peak base shear of frame B (kips).
VbBSRSS=SRSS(VbB);
%% Modal combination with the Complete Quadratic Combination (CQC) method.
% Peak displacement in y direction (inch).
uyCQC=CQC(U(2,:)',omega,ksi);
%%
% Peak torsional rotation (rad).
uthCQC=CQC(U(3,:)',omega,ksi);
%%
% Peak base shear of the whole structure (kips).
VbCQC=CQC(f(2,:)',omega,ksi);
%%
% Peak base torsional moment of the whole structure (kips*ft).
TbCQC=CQC(f(3,:)',omega,ksi);
%%
% Peak base shear of frame A (kips).
VbACQC=CQC(VbA,omega,ksi);
%%
% Peak base shear of frame B (kips).
VbBCQC=CQC(VbB,omega,ksi);
%% Assemble values of peak response in a table.
% Assemble values of peak response in a cell.
C{1,2}='uy (in)';
C{1,3}='(b/2)uth (in)';
C{1,4}='Vb (kips)';
C{1,5}='Tb (kip-ft)';
C{1,6}='VbA (kips)';
C{1,7}='VbB (kips)';
C{2,1}='ABSSUM';
C{2,2}=uyAbsSum;
C{2,3}=(b/2)*uthAbsSum*12;
C{2,4}=VbAbsSum;
C{2,5}=TbAbsSum;
C{2,6}=VbAAbsSum;
C{2,7}=VbBAbsSum;
C{3,1}='SRSS';
C{3,2}=uySRSS;
C{3,3}=(b/2)*uthSRSS*12;
C{3,4}=VbSRSS;
C{3,5}=TbSRSS;
C{3,6}=VbASRSS;
C{3,7}=VbBSRSS;
C{4,1}='CQC';
C{4,2}=uyCQC;
C{4,3}=(b/2)*uthCQC*12;
C{4,4}=VbCQC;
C{4,5}=TbCQC;
C{4,6}=VbACQC;
C{4,7}=VbBCQC;
C
%% Verification of table E13.12b of Chopra (Dynamics of Structures, 2012)
% Compare the peak responses with those shown in table E13.12b of Chopra (2012).
%
% <<ChopraE1312b.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%



