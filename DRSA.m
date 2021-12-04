function [U,V,A,f,omega,Eigvec] = DRSA(K,M,r,PSa,Sd,Sv,Sa,eigInd)
%
% Peak modal structural response
%
% function [U,V,A,F,OMEGA,EIGVEC] = DRSA(K,M,R,PSA,SD,SV,SA,EIGIND)
%
% Description
%     Determine the peak structural response using the linear elastic
%     response spectrum of the ground motion, without carrying a response
%     history analysis.
%
% Input parameters
%     K [double(:NDOFS x :NDOFS)]: Stiffness matrix of the structure to
%         be analysed. #ndofs# is the number of degrees of freedom of the
%         structure
%     M [double(:NDOFS x :NDOFS)]: Mass matrix of the structure to be
%         analysed
%     R [double(:NDOFS x 1)]: Influence vector. It determines the spatial
%         distribution of the effective earthquake forces
%     PSA [double(:inf x 2)]: Pseudo spectral acceleration. The first
%         column contains the period values and the second column contains
%         the pseudospectral acceleration values.
%     SD [double(:inf x 2)]: Spectral displacement. The first column
%         contains the period values and the second column contains the
%         spectral displacement values.
%     SV [double(:inf x 2)]: Spectral velocity. The first column contains
%         the period values and the second column contains the spectral
%         velocity values.
%     SA [double(:inf x 2)]: Spectral acceleration. The first column
%         contains the period values and the second column contains the
%         spectral acceleration values.
%     EIGIND [double(:inf x 1)]: Eigenmode indicator. Only the eigenmode
%         numbers that are contained in EIGIND are taken into account for
%         the calculation of the peak modal responses.
%
% Output parameters
%     U [double(1:NDOFS x 1:NEIG)]: Peak modal displacements.
%     V [double(1:NDOFS x 1:NEIG)]: Peak modal velocities
%     A [double(1:NDOFS x 1:NEIG)]: Peak modal accelerations
%     F [double(1:NDOFS x 1:NEIG)]: Peak modal equivalent static forces
%     OMEGA [double(1:NEIG x 1)]: Cyclic eigenfrequencies
%     EIGVEC [double(1:NDOFS x 1:NEIG)]: Eigenvectors (stored as columns)
%         sorted in ascending order of eigenfrequency
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

ndofs=size(K,1);
neig=numel(eigInd);
% Calculate eigenvalues, eigenvectors and their number
[Eigvec,Eigval]=eig(K,M);
Eigvec=Eigvec(:,eigInd);
% Take the eigenvalues in column vector and sort in ascending order of
% eigenfrequency:
D1=diag(Eigval,0);
D1=D1(eigInd);
% Take the eigenvectors (stored as columns) and sort in ascending order of
% eigenfrequency:
%[D2,IX]=sort(D1);
%EigvecSort=Eigvec(:,IX);
% Generalized masses Mn for all eigenmodes from eq.(13.1.5) of Chopra
% (2012).
Mn=diag(Eigvec'*M*Eigvec);
% Ln coefficients from eq.(13.1.5) of Chopra (2012).
Ln=Eigvec'*M*r;
% Gamman coefficients from eq.(13.1.5) of Chopra (2012).
Gamman=Ln./Mn;
% Eigenperiods of the structure
omega=D1.^0.5;
Ti=2*pi./omega;

% Peak modal displacements
if ~isempty(Sd)
    Sdi = interp1(Sd(:,1),Sd(:,2),Ti,'linear','extrap');
    U=zeros(ndofs,neig);
    for i=1:neig
        U(:,i)=Gamman(i)*Eigvec(:,i)*Sdi(i);
    end
else
    U=[];
end
% Peak modal velocities
if ~isempty(Sv)
    Svi = interp1(Sv(:,1),Sv(:,2),Ti,'linear','extrap');
    V=zeros(ndofs,neig);
    for i=1:neig
        V(:,i)=Gamman(i)*Eigvec(:,i)*Svi(i);
    end
else
    V=[];
end
% Peak modal accelerations
if ~isempty(Sa)
    Sai = interp1(Sa(:,1),Sa(:,2),Ti,'linear','extrap');
    A=zeros(ndofs,neig);
    for i=1:neig
        A(:,i)=Gamman(i)*Eigvec(:,i)*Sai(i);
    end
else
    A=[];
end
% Peak modal equivalent static forces
if ~isempty(PSa)
    PSai = interp1(PSa(:,1),PSa(:,2),Ti,'linear','extrap');
    f=zeros(ndofs,neig);
    for i=1:neig
        f(:,i)=Gamman(i)*(M*Eigvec(:,i))*PSai(i);
    end
else
    f=[];
end


end

