% REMRND        round to rem.
%
% Y = REMRND( X, REM, MODE )
%
% ex: remrnd( 7.2, 2 ) returns 8.
%     remrnd( 7.2, 2, 'ceil' ) returns 8
%     remrnd( 7.2, 2, 'floor' ) returns 6

% 24-jan-04 ES

% 28-feb-12 tol added for numerical stability

function y = remrnd( x, rem, mode, tol )

if nargin < 3 || isempty( mode ), mode = 'round'; end
if nargin < 4 || isempty( tol ), tol = sqrt( eps ); end

switch mode
    case 'ceil'
        y = ceil( ( x - tol ) ./ rem ) .* rem;
    case 'floor'
        y = floor( ( x + tol ) ./ rem ) .* rem;
    otherwise
        y = round( x ./ rem ) .* rem;
end

return