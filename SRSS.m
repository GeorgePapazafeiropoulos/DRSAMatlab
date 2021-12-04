function y = SRSS(x)
%
% Modal combination with the Square Root of Sum of Squares (SRSS) method
%
% function Y = SRSS(X)
%
% Input parameters
%     Y [double(:inf x 1)]: A column vector containing the peak response
%         for each eigenmode of a structure.
%
% Output parameters
%     X [double(1 x 1)]: Combination of the peak responses of the various
%         eigenmodes according to the Square Root of the Sum of Squares
%         (SRSS) rule.
%
% Example
%     x=rand(5,1)-0.5;
%     y=SRSS(x)
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

y=sqrt(sum(x.^2));

end

