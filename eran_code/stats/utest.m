% UTEST         Mann-Whitney U-test.
%
% call          [ P U METHOD ] = UTEST( X1, X2, TEST )
%
% gets          X1, X2      samples; may be unequal size
%               TEST        test type:
%                               {0}     two sided, H1: X1 != X2
%                                1      one sided, H1: X1 >  X2
%                               -1      one sided, H1: X1 <  X2
%
% returns       P           p-value (of H0: X1 = X2)
%               U           the statistic
%               METHOD      method used:
%                               N1+N2 small, no ties        table lookup
%                               N1+N2 small, ties           exact binomial
%                               N1+N2 large                 Gaussian approximation
%
% calls         RANKCOLS, UTABLE
%
% see also      MAKE_UTABLE
%
% note          the numerical results obtained by this function should
%               be identical to a two-sided Wilcoxon's rank-sum test; this
%               routine provides 2 advantages: (1) one sided testing (2)
%               table lookup in small sample cases

% reference: Sokal & Rohlf 2001, p.427-431

% 18-jul-04 ES

% revisions
% 03-may-07 support of NaNs
% 03-aug-07 make sure column vectors
% 22-apr-13 support single 2-col matrix input
% 12-sep-19 cleaned up

function [ p, U, method ] = utest( x1, x2, test )

% handle arguments

nargs = nargin;
if nargs < 2
    if size( x1, 2 ) == 2
        x2 = x1( :, 2 );
        x1 = x1( :, 1 );
    else
        error( '2 arguments' )
    end
end
if nargs < 3 || isempty( test ), test = 0; end

% handle special cases

p = NaN;
U = NaN;
method = 'none'; 
if isequal( x1, x2 ), p = 1; return, end
if isempty( x1 ) || isempty( x2 ), p = 1; return, end

% sample sizes
x1 = x1( ~isnan( x1 ) );
x2 = x2( ~isnan( x2 ) );
x1 = x1( : );
x2 = x2( : );
n1 = length( x1 ); 
n2 = length( x2 );
n = n1 + n2;
if n1 == 0 || n2 == 0
    return
end
    
% Mann-Whitney U statistic

[ ranks, T ] = rankcols( [ x1; x2 ] );
R = sum( ranks( 1 : n1 ) );
C = n1 * n2 + n1 * ( n1 + 1 ) / 2 - R;              % p.429
U = max( C, n1 * n2 - C );

% significance

expcR = n1 * ( n + 1 ) / 2;
if ( n1 + n2 ) < 20
    if T                                            % exact binomial
        f = sum( nchoosek( ranks, n1 ), 2 );
        if R < expcR
            p = sum( f <= R ) / length( f );
        else 
            p = sum( f >= R ) / length( f );
        end
        method = 'exact';
    else                                            % use a table based on no-ties (perfect for that case only)
        eval( sprintf( 'load( ''utable.mat'', ''ptab_%d'' )', n ) )
        if R < expcR
            eval( sprintf( 'p = ptab_%d( ceil( R ), n1 );', n ) )
        else
            eval( sprintf( 'p = 1 - ptab_%d( ceil( R ) - 1, n1 );', n ) )
        end
        method = 'lookup';
    end
else                                                % approximate gaussian
    nom = U - n1 * n2 / 2;
    if T                                            % correction for ties, p.430
        den = n1 * n2 / 12 / n / ( n - 1 ) * ( n .^ 3 - n - T );
    else
        den = n1 * n2 / 12 * ( n + 1 );
    end
    ts = nom / sqrt( den );
    p = normcdf( -abs( ts ), 0, 1 );
    method = 'gaussian';
end

% translate p-value according to test type

if ~test
    p = 2 * p;
elseif ( test == 1 && R < expcR ) || ( test == -1 && R > expcR )
    p = 1 - p;
end

return

% to test (and compare with matlab's ranksum)

n1 = 8; n2 = 6; 
x1 = round( rand( 10, 1 ) * 100 ) + 30; x2 = round( rand( 6, 1 ) * 100 ) + 5;   % generate random samples
t0 = clock; p_rs = ranksum( x1, x2 ); t( 1 ) = etime( clock, t0 ); 
t0 = clock; [ p_u U, method ] = utest( x1, x2 ); t( 2 ) = etime( clock, t0 ); 
fprintf( 1, '%s\nn: [ %d %d ]; medians: [ %0.3g %0.3g ]\n', repmat( '*', 1, 80 ), n1, n2, median( x1 ), median( x2 ) );
fprintf( 1, 'rs: %0.4g, t = %0.3g\nut: %0.4g, t = %0.3g, method: %s\n%s\n', p_rs, t( 1 ), p_u, t( 2 ), method, repmat( '*', 1, 80 ) );
