% BOUNDS            lower and upper bounds
%
% call              [ B MSG ] = BOUNDS( X, ALPHA )
%
% gets              X       data. if matrix, works on columns
%                   ALPHA   1-confidence 
%                           if paired, uses non-symmetric limits. 
%                           if a column vector (or a matrix), evaluates
%                           multiple bounds for a vector X (if X is a
%                           matrix, only first ALPHA is evaluated)
%
% returns           B       [ LB; UB ]
%                   MSG     warnings (nans, insufficient data)
%
% calls             nothing

% 27-mar-03 ES

% revisions
% 22-jun-03 edges warning
% 04-sep-03 remove nans before computing
% 16-dec-03 msg returned
% 02-may-04 edge conditions improved; non-symmetric and multiple alphas supported
% 18-jul-04 matrices infested with NaN handled properly
% 15-dec-19 cleaned up

function [ b, msg ] = bounds( x, alpha )

msg                     = '';
if nargin < 2, alpha    = 0.01; end

% handle different input formats
sx                      = size( x );
sa                      = size( alpha );
if any( sx == 1 )                               % single vector input
    x                   = x( : );
    mat                 = 0;
else                                            % matrix input
    mat                 = 1;
end
if sa( 2 ) == 1                                 % symmetric bounds
    alpha               = [ alpha / 2, 1 - alpha / 2 ];
end
if mat && sa( 1 ) ~= 1                           % matrix input AND multiple alphas
    alpha               = alpha( 1, : );
    sa( 1 )             = 1;
end

% get the actual bounds
xs                      = sort( x );
lx                      = sum( ~isnan( x ) );
si                      = max( floor( alpha( :, 1 ) * lx ), 1 );
ei                      = min( ceil( alpha( :, 2 ) * lx ), ones( sa( 1 ), 1 ) * lx );

if mat
    si                  = si + ( 0 : sx( 1 ) : sx( 1 ) * ( sx( 2 ) - 1 ) );
    ei                  = ei + ( 0 : sx( 1 ) : sx( 1 ) * ( sx( 2 ) - 1 ) );
    b                   = [ xs( si ); xs( ei ) ];
else
    b                   = [ xs( si )'; xs( ei )' ];
end

% report ill conditions
if any( si == 1 ) || any( ei == lx )
    msg                 = sprintf( '%sinsufficient data; ', msg );
end
if any( lx ~= sx( 1 ) )
    msg                 = sprintf( '%sNaNs; ', msg );
end

return

% EOF
