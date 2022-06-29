% effective_independence    number of independent measurements from data
% 
% call                      [ m, w ] = effective_independence( x )
%
% gets                      x       matrix with data in columns
%
% returns                   m       number of independent measurements
%                           w       full-width at half-height of each data set
%
% calls                     calc_fwhh, my_xcorr (ssp)

% 30-dec-19 ES

function [ m, w ] = effective_independence( x )

m               = NaN;
if nargin < 1 || isempty( x )
    return
end
ac              = my_xcorr( x ); 
w               = calc_fwhh( ac ); 
n               = size( x, 1 );
m               = n ./ w;

return

% EOF