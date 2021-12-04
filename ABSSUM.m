function y = ABSSUM(x)
%
% Modal combination with the ABSolute SUM (ABSSUM) method
%
% function Y = ABSSUM(X)
%
% Input parameters
%     X [double(:inf x 1)]: A column vector containing the peak response
%         for each eigenmode of a structure.
%
% Output parameters
%     Y [double(1 x 1)]: Combination of the peak responses of the various
%         eigenmodes according to the sum of the absolute (ABSSUM) values
%         rule.
%
% Example
%     x=rand(5,1)-0.5;
%     y=ABSSUM(x)
%
%__________________________________________________________________________
% Copyright (c) 2015-2021
%     George Papazafeiropoulos
%     Major, Infrastructure Engineer, Hellenic Air Force
%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%     Email: gpapazafeiropoulos@yahoo.gr
% _________________________________________________________________________

y=sum(abs(x));

end

