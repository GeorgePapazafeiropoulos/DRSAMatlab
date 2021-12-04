%% 2-DOF shear buiding subjected to motion with exponential response spectrum

%% Statement of the problem
% * A 2-degree-of-freedom system is subjected to horizontal earthquake
% excitation (figure below) with its response spectra as given below. Take
% the flexural rigidity, EI = 10^6 Nm^2 and length, L = 2m. Each lumped
% mass is 100 kg. Determine the maximum top mass floor displacement and
% base shear. Take 2% damping in each mode of vibration.
%
% <<Shear_Building_2DOF_ExpSpec_fig1.png>>
%
% <<Shear_Building_2DOF_ExpSpec_fig2.png>>
%
%% Initialization of structural input data
% Set the storey height of the structure in m.
L=2;
%%
% Set the number of eigenmodes of the structure, which is equal to the
% number of its storeys.
neig=2;
%%
% Set the flexural rigidity of each floor in Nm^2
EI=1e6;
%%
% Set the lumped mass at each floor in kg.
m=100;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=6*EI/(7*L^3)*[16,-5;-5,2];
%%
% Calculate the mass matrix of the structure.
M=m*eye(neig);
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=ones(2,1);
%% Load earthquake response spectrum data
% Earthquake pseudoacceleration response spectrum of the El Centro
% earthquake (El Centro, 1940, El Centro Terminal Substation Building)
T=linspace(0.01,5,100)';
PSa=[T,9.81*(0.1+T).*exp(-T)];
%%
% Earthquake displacement response spectrum of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
Sd=[T,PSa(:,2)./((2*pi./T).^2)];
%%
% Set the critical damping ratio of the response spectra
% ($$\mathrm{\xi}=0.05$)
ksi=0.02;
%% Dynamic Response Spectrum Analysis (DRSA)
% Consider all eigenmodes of the building
eigInd=[1;2];
%%
% Perform DRSA analysis
[U,~,~,f,omega,Eigvec] = DRSA(K,M,r,PSa,Sd,[],[],eigInd);
%%
% Plot the natural modes of vibration of the uniform five-story shear
% building
FigHandle=figure('Name','Natural Modes','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;Eigvec(:,i)],(0:L:L*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 0 1],'markeredgecolor','k')
    grid on
    xlabel('Displacement','FontSize',13);
    ylabel('Height (m)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Plot the peak modal displacement response.
FigHandle=figure('Name','Displacements','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;U(:,i)],(0:L:L*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 1 0],'markeredgecolor','k')
    xlim([-max(abs(U(:,i))) max(abs(U(:,i)))])
    grid on
    xlabel('Displacement (m)','FontSize',13);
    ylabel('Height (m)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Plot the peak modal equivalent static force response.
FigHandle=figure('Name','Equivalent static forces','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;f(:,i)],(0:L:L*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[1 0 0],'markeredgecolor','k')
    xlim([-max(abs(f(:,i))) max(abs(f(:,i)))])
    grid on
    xlabel('Static force (N)','FontSize',13);
    ylabel('Height (m)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Calculate the peak modal base shear in N.
Vb=zeros(neig,1);
for i=1:neig
    Vb(i)=sum(f(:,i));
end
Vb
%%
% Calculate the peak modal base overturning moment in Nm.
Mb=zeros(neig,1);
for i=1:neig
    Mb(i)=sum(f(:,i).*(L:L:2*L)');
end
Mb
%% Modal combination with the ABSolute SUM (ABSSUM) method.
% Peak base shear.
VbAbsSum=ABSSUM(Vb);
%%
% Peak top-story shear.
V2AbsSum=ABSSUM(f(2,:)');
%%
% Peak base overturning moment.
MbAbsSum=ABSSUM(Mb);
%%
% Peak top-story displacement.
u2AbsSum=ABSSUM(U(2,:)');
%% Modal combination with the Square Root of Sum of Squares (SRSS) method.
% Peak base shear.
VbSRSS=SRSS(Vb);
%%
% Peak top-story shear.
V2SRSS=SRSS(f(2,:)');
%%
% Peak base overturning moment.
MbSRSS=SRSS(Mb);
%%
% Peak top-story displacement.
u2SRSS=SRSS(U(2,:)');
%% Modal combination with the Complete Quadratic Combination (CQC) method.
% Peak base shear.
VbCQC=CQC(Vb,omega,ksi);
%%
% Peak top-story shear.
V2CQC=CQC(f(2,:)',omega,ksi);
%%
% Peak base overturning moment.
MbCQC=CQC(Mb,omega,ksi);
%%
% Peak top-story displacement.
u2CQC=CQC(U(2,:)',omega,ksi);
%% Assemble values of peak response in a table.
% Assemble values of peak response in a cell.
C{1,2}='Vb (N)';
C{1,3}='V2 (N)';
C{1,4}='Mb (Nm)';
C{1,5}='u2 (m)';
C{2,1}='ABSSUM';
C{2,2}=VbAbsSum;
C{2,3}=V2AbsSum;
C{2,4}=MbAbsSum;
C{2,5}=u2AbsSum;
C{3,1}='SRSS';
C{3,2}=VbSRSS;
C{3,3}=V2SRSS;
C{3,4}=MbSRSS;
C{3,5}=u2SRSS;
C{4,1}='CQC';
C{4,2}=VbCQC;
C{4,3}=V2CQC;
C{4,4}=MbCQC;
C{4,5}=u2CQC;
C
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%


