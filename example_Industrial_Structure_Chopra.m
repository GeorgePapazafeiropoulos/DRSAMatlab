%% Industrial structure subjected to the El Centro earthquake motion (Chopra, 2012)

%% Statement of the problem
% An industrial structure is modeled as 2-DOF system as shown in the
% following figure. Determine the lateral displacement, base shear and base
% moment of the structure due to El-Centro, 1940 earthquake ground motion
% using the response spectrum method. Take $$EI=80000 Nm^2$, $$L=2m$,
% $$m_1=100kg$ and $$m_2=200kg$. The damping ratio for all natural modes
% is $$\mathrm{\xi} = 0.02$.
%
% <<Industrial.png>>
%
%% Initialization of structural input data
% Set the length of the structure members in m.
L=2;
%%
% Set the flexural stiffness in Nm^2.
EI=80000;
%%
% Set the lumped masses in kg.
m1=100;
m2=200;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=6*EI/(7*L^3)*[8,-3;-3,2];
%%
% Calculate the mass matrix of the structure in kg.
M=[m1+m2,0;0,m2];
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied only at dof No 1 of the structure.
r=[1;0];
%% Load earthquake response spectrum data
% Earthquake displacement response spectrum of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
Sd=load('Sd2.txt');
%%
% Earthquake pseudoacceleration response spectrum of the El Centro
% earthquake (El Centro, 1940, El Centro Terminal Substation Building)
PSa=load('PSa2.txt');
%%
% Set the critical damping ratio of the response spectra
% ($$\mathrm{\xi}=0.02$)
ksi=0.02;
%% Dynamic Response Spectrum Analysis (DRSA)
% Consider all eigenmodes of the building
eigInd=[1;2];
%%
% Perform DRSA analysis
[U,~,~,f,omega,Eigvec] = DRSA(K,M,r,PSa,Sd,[],[],eigInd);
%%
% Set the number of eigenmodes of the structure.
neig=numel(omega);
%%
% Plot the natural modes of vibration of the industrial structure
FigHandle=figure('Name','Natural Modes','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 1000, 500]);
for i=1:neig
    subplot(neig,1,i)
    plot([0;10*Eigvec(1,i);L],[0;L;10*Eigvec(2,i)+L],'LineWidth',2.,...
        'Marker','.','MarkerSize',20,'Color',[0 0 1],'markeredgecolor','k')
    grid on
    xlabel('Displacement','FontSize',13);
    ylabel('Height','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Calculate the peak modal base shear (N).
Vb=zeros(neig,1);
for i=1:neig
    Vb(i)=f(1,i);
end
Vb
%%
% Calculate the peak modal base overturning moment (Nm).
Mb=zeros(neig,1);
for i=1:neig
    Mb(i)=L*sum(f(:,i));
end
Mb
%% Modal combination with the ABSolute SUM (ABSSUM) method.
% Calculation of peak base shear.
VbAbsSum=ABSSUM(Vb);
%%
% Calculation of peak base overturning moment.
MbAbsSum=ABSSUM(Mb);
%%
% Calculation of peak lateral displacement.
u1AbsSum=ABSSUM(U(1,:)');
%% Modal combination with the Square Root of Sum of Squares (SRSS) method.
% Calculation of peak base shear.
VbSRSS=SRSS(Vb);
%%
% Calculation of peak base overturning moment.
MbSRSS=SRSS(Mb);
%%
% Calculation of peak lateral displacement.
u1SRSS=SRSS(U(1,:)');
%% Modal combination with the Complete Quadratic Combination (CQC) method.
% Calculation of peak base shear.
VbCQC=CQC(Vb,omega,ksi);
%%
% Calculation of peak base overturning moment.
MbCQC=CQC(Mb,omega,ksi);
%%
% Calculation of peak lateral displacement.
u1CQC=CQC(U(1,:)',omega,ksi);
%% Assemble values of peak response in a table.
% 
C{1,2}='Lateral Displ. (m)';
C{1,3}='Base Shear (N)';
C{1,4}='Base Moment (Nm)';
C{2,1}='ABSSUM';
C{2,2}=u1AbsSum;
C{2,3}=VbAbsSum;
C{2,4}=MbAbsSum;
C{3,1}='SRSS';
C{3,2}=u1SRSS;
C{3,3}=VbSRSS;
C{3,4}=MbSRSS;
C{4,1}='CQC';
C{4,2}=u1CQC;
C{4,3}=VbCQC;
C{4,4}=MbCQC;
C


%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%


