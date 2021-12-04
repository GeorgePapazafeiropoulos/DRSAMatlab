%% 2-DOF shear buiding analyzed according to a design response spectrum

%% Statement of the problem
% * For the design spectrum below calculate the maximum displacements and
% interstorey drifts of the two floors and the base shear. Take behavior
% factor equal to 3.
%
% <<Shear_Building_2DOF_DesSpec_fig1.png>>
%
% <<Shear_Building_2DOF_DesSpec_fig2.png>>
%
%% Initialization of structural input data
% Set the storey stiffnesses of the structure in N/m.
k1=200e6;
k2=150e6;
%%
% Set the storey lumped masses of the structure in kg.
m1=100e3;
m2=80e3;
%%
% Set the number of eigenmodes of the structure, which is equal to the
% number of its storeys.
neig=2;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in N/m.
K=[k2,-k2;-k2,k2+k1];
%%
% Calculate the mass matrix of the structure.
M=[m2,0;0,m1];
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=[1;1];
%% Load earthquake response spectrum data
% Pseudoacceleration design response spectrum
T1=linspace(0.01,0.24,30)';
T2=linspace(0.25,5,100)';
PSa=[[T1;T2],[0.3*9.81*ones(size(T1));0.3*9.81*(0.25./T2).^(2/3)]];
%%
% Displacement design response spectrum
Sd=[PSa(:,1),PSa(:,2).*(PSa(:,1).^2)/4/pi^2];
%%
% Set the behavior factor
q=3;
%% Dynamic Response Spectrum Analysis (DRSA)
% Consider all eigenmodes of the building
eigInd=[1;2];
%%
% Perform DRSA analysis
[U,~,~,f,omega,Eigvec] = DRSA(K,M,r,PSa,Sd,[],[],eigInd);
%%
% Cyclic eigenfrequencies of the structure
omega
%%
% Normalize the eigenmodes of the structure
phi1=Eigvec(:,1)/Eigvec(1,1)
phi2=Eigvec(:,2)/Eigvec(1,2)
%%
% Design forces
f_design=f/q
%%
% Plot the natural modes of vibration of the uniform five-story shear
% building
FigHandle=figure('Name','Natural Modes','NumberTitle','off');
set(FigHandle,'Position',[50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;Eigvec(:,i)],(0:1:neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 0 1],'markeredgecolor','k')
    grid on
    xlabel('Displacement','FontSize',13);
    ylabel('Height','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Plot the peak modal displacement response.
FigHandle=figure('Name','Displacements','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;U(:,i)],(0:1:neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 1 0],'markeredgecolor','k')
    xlim([-max(abs(U(:,i))) max(abs(U(:,i)))])
    grid on
    xlabel('Displacement (m)','FontSize',13);
    ylabel('Height','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Plot the peak modal equivalent static force response.
FigHandle=figure('Name','Equivalent static forces','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;f_design(:,i)],(0:1:neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[1 0 0],'markeredgecolor','k')
    xlim([-max(abs(f_design(:,i))) max(abs(f_design(:,i)))])
    grid on
    xlabel('Static force (N)','FontSize',13);
    ylabel('Height','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Calculate the design modal base shear in N.
Vb=zeros(neig,1);
for i=1:neig
    Vb(i)=sum(f_design(:,i));
end
Vb
%% Modal combination with the Square Root of Sum of Squares (SRSS) method.
% Design base shear.
VbSRSS=SRSS(Vb)
%%
% Design top-story shear.
V2SRSS=SRSS(f_design(2,:)')
%%
% Design top-story displacement.
u2SRSS=SRSS(U(2,:)')
%%
% Design interstorey drift between 1st and 2nd storey.
i21SRSS=SRSS(U(2,:)'-U(1,:)')
%%
% Design interstorey drift between 1st storey and ground.
u1SRSS=SRSS(U(1,:)')

%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%


