% calc_fwhh           determine full-width at half-height (direct method)
%
% [ width lag val ] = calc_fwhh( x )
%
% if x is a matrix, works on columns
% 
% width: FWHH [samples]
% lag: of extremum [samples]
% val: of the extremum
% 
% NOTE: assumes a zero reference

% 23-mar-14 ES

function [ width lag val ] = calc_fwhh( x )

[ m n ] = size( x );
width = NaN * ones( 1, n );
lag = width;
val = width;
idx = NaN * ones( m, n );

if all( isnan( x ) )
    return
end
    
ref = zeros( 1, n );
midx = max( abs( x ), [], 1 ) == max( x, [], 1 );
if sum( midx ) > 0
    [ val( midx ) lag( midx ) ] = max( x( :, midx ) );
    halfAmp = ( val( midx ) + ref( midx ) ) / 2;
    idx( :, midx ) = cumsum( bsxfun( @lt, x( :, midx ), halfAmp ) );
end
midx = max( abs( x ), [], 1 ) == -min( x, [], 1 );
if sum( midx ) > 0
    [ val( midx ) lag( midx ) ] = min( x( :, midx ) );
    halfAmp = ( val( midx ) + ref( midx ) ) / 2;
    idx( :, midx ) = cumsum( bsxfun( @gt, x( :, midx ), halfAmp ) );
end

% eidx = idx( lag + m * ( 0 : ( n - 1 ) ) );
% width = sum( bsxfun( @eq, idx, eidx ) ) - 1;

nnans = ~isnan( lag );
idxvalid = idx( :, nnans );
eidx = idxvalid( lag( nnans ) + m * ( 0 : ( sum( nnans ) - 1 ) ) );
width( nnans ) = sum( bsxfun( @eq, idxvalid, eidx ) ) - 1;

return

% EOF
% 9    17    21    24    29    32    34    37    41    47