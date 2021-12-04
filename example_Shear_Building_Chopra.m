%% Shear buiding subjected to the El Centro earthquake motion (Chopra, 2012)

%% Statement of the problem
% * *Chopra (2012), Section 12.8:* Consider the structure of Fig. 12.8.1: a
% uniform five-story shear building (i.e., flexurally rigid floor beams and
% slabs) with lumped mass m at each floor, and same story stiffness k for
% all stories.
% * *Chopra (2012), Section 13.2.6:* The structure is subjected to the El
% Centro ground motion (Chopra (2012), Fig. 6.1.4). The lumped mass $$m_j =
% m = 100 kips/g$ at each floor, the lateral stiffness of each story is
% $$k_j = k = 31.54 kips/in.$, and the height of each story is 12 ft. The
% damping ratio for all natural modes is $$\mathrm{\zeta_n} = 0.05$.
%
% <<Chopra1281.png>>
%
%% Initialization of structural input data
% Set the storey height of the structure in ft.
h=12;
%%
% Set the number of eigenmodes of the structure, which is equal to the
% number of its storeys.
neig=5;
%%
% Set the lateral stiffness of each storey in kips/inch.
k=31.54;
%%
% Set the lumped mass at each floor (g=386.4 inch/sec^2).
m=100/9.81*0.0254;
%% Calculation of structural properties
% Calculate the stiffness matrix of the structure in kips/inch.
K=k*(diag([2*ones(neig-1,1);1])+diag(-ones(neig-1,1),1)+diag(-ones(neig-1,1),-1));
%%
% Calculate the mass matrix of the structure.
M=m*eye(neig);
%%
% Set the spatial distribution of the effective earthquake forces.
% Earthquake forces are applied at all dofs of the structure.
r=ones(5,1);
%% Load earthquake response spectrum data
% Earthquake displacement response spectrum of the El Centro earthquake (El
% Centro, 1940, El Centro Terminal Substation Building)
Sd=load('Sd1.txt');
%%
% Earthquake pseudoacceleration response spectrum of the El Centro
% earthquake (El Centro, 1940, El Centro Terminal Substation Building)
PSa=load('PSa1.txt');
%%
% Set the critical damping ratio of the response spectra
% ($$\mathrm{\xi}=0.05$)
ksi=0.05;
%% Dynamic Response Spectrum Analysis (DRSA)
% Consider all eigenmodes of the building
eigInd=[1;2;3;4;5];
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
    plot([0;Eigvec(:,i)],(0:h:h*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 0 1],'markeredgecolor','k')
    grid on
    xlabel('Displacement','FontSize',13);
    ylabel('Height (ft)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Compare the eigenmodes with those shown in Figure 12.8.2. in Chopra
% (2012).
%
% <<Chopra1282.png>>
%
%%
% Plot the peak modal displacement response.
FigHandle=figure('Name','Displacements','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;U(:,i)],(0:h:h*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[0 1 0],'markeredgecolor','k')
    xlim([-max(abs(U(:,i))) max(abs(U(:,i)))])
    grid on
    xlabel('Displacement (in)','FontSize',13);
    ylabel('Height (ft)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Plot the peak modal equivalent static force response.
FigHandle=figure('Name','Equivalent static forces','NumberTitle','off');
set(FigHandle, 'Position', [50, 50, 1000, 500]);
for i=1:neig
    subplot(1,neig,i)
    plot([0;f(:,i)],(0:h:h*neig)','LineWidth',2.,'Marker','.',...
        'MarkerSize',20,'Color',[1 0 0],'markeredgecolor','k')
    xlim([-max(abs(f(:,i))) max(abs(f(:,i)))])
    grid on
    xlabel('Static force (kips)','FontSize',13);
    ylabel('Height (ft)','FontSize',13);
    title(['Mode ',num2str(i)],'FontSize',13)
end
%%
% Calculate the peak modal base shear in kips.
Vb=zeros(neig,1);
for i=1:neig
    Vb(i)=sum(f(:,i));
end
Vb
%%
% Calculate the peak modal base overturning moment in kips-ft.
Mb=zeros(neig,1);
for i=1:neig
    Mb(i)=sum(f(:,i).*(h:h:5*h)');
end
Mb
%% Verification of figure 13.8.3 of Chopra (Dynamics of Structures, 2012)
% Compare the peak modal responses with those shown in figure 13.8.3 of
% Chopra (2012).
%
% <<Chopra1383.png>>
%
%% Modal combination with the ABSolute SUM (ABSSUM) method.
% Peak base shear.
VbAbsSum=ABSSUM(Vb);
%%
% Peak top-story shear.
V5AbsSum=ABSSUM(f(5,:)');
%%
% Peak base overturning moment.
MbAbsSum=ABSSUM(Mb);
%%
% Peak top-story displacement.
u5AbsSum=ABSSUM(U(5,:)');
%% Modal combination with the Square Root of Sum of Squares (SRSS) method.
% Peak base shear.
VbSRSS=SRSS(Vb);
%%
% Peak top-story shear.
V5SRSS=SRSS(f(5,:)');
%%
% Peak base overturning moment.
MbSRSS=SRSS(Mb);
%%
% Peak top-story displacement.
u5SRSS=SRSS(U(5,:)');
%% Modal combination with the Complete Quadratic Combination (CQC) method.
% Peak base shear.
VbCQC=CQC(Vb,omega,ksi);
%%
% Peak top-story shear.
V5CQC=CQC(f(5,:)',omega,ksi);
%%
% Peak base overturning moment.
MbCQC=CQC(Mb,omega,ksi);
%%
% Peak top-story displacement.
u5CQC=CQC(U(5,:)',omega,ksi);
%% Assemble values of peak response in a table.
% Assemble values of peak response in a cell.
C{1,2}='Vb (kips)';
C{1,3}='V5 (kips)';
C{1,4}='Mb (kip-ft)';
C{1,5}='u5 (in)';
C{2,1}='ABSSUM';
C{2,2}=VbAbsSum;
C{2,3}=V5AbsSum;
C{2,4}=MbAbsSum;
C{2,5}=u5AbsSum;
C{3,1}='SRSS';
C{3,2}=VbSRSS;
C{3,3}=V5SRSS;
C{3,4}=MbSRSS;
C{3,5}=u5SRSS;
C{4,1}='CQC';
C{4,2}=VbCQC;
C{4,3}=V5CQC;
C{4,4}=MbCQC;
C{4,5}=u5CQC;
C
%% Verification of table 13.8.5 of Chopra (Dynamics of Structures, 2012)
% Compare the peak responses with those shown in table 13.8.5 of Chopra
% (2012).
%
% <<Chopra1385.png>>
%
%% Copyright
%
% Copyright (c) 2015-2021 by George Papazafeiropoulos
%
% * Major, Infrastructure Engineer, Hellenic Air Force
% * Civil Engineer, M.Sc., Ph.D. candidate, NTUA
% * Email: gpapazafeiropoulos@yahoo.gr
%


