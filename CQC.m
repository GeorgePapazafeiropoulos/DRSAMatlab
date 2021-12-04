function y = CQC(x,omega,ksi)
%
% Modal combination with the Complete Quadratic Combination (CQC) method
%
% function Y = CQC(X,OMEGA,KSI)
%
% Input parameters
%     X [double(:inf x 1)]: A column vector containing the peak response
%         for each eigenmode of a structure.
%     OMEGA [double(:inf x 1)]: Cyclic eigenfrequencies of the eigenmodes
%         of the structure. SIZE(OMEGA,1) must be equal to SIZE(X,1).
%     KSI [double(1 x 1)]: Critical damping ratio
%
% Output parameters
%     Y [double(1 x 1)]: Combination of the peak responses of the various
%         eigenmodes according to the Complete Quadratic Combination (CQC)
%         rule.
%
% Example
%     x=rand(5,1)-0.5;
%     omega=rand(5,1);
%     ksi=0.05;
%     y = CQC(x,omega,ksi)
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

% Ratios of cyclic frequencies.
neig=numel(omega);
bij=omega(:,ones(1,neig))./omega(:,ones(1,neig))';
% Correlation coefficients:
aij=8.*ksi.^2.*(1+bij).*bij.^(3/2)./((1-bij.^2).^2+4.*ksi.^2.*bij.*(1+bij).^2);
% Matrix with double summation terms:
x=x';
xij=x(ones(neig,1),:).*x(ones(neig,1),:)'.*aij;
% Sum all terms and take square root
y=sqrt(sum(xij(:)));

end

