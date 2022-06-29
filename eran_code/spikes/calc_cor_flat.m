% CALC_COR_FLAT         calculate asymptotically flat CCH
%
% call                  CCH = CALC_COR_FLAT( ST1, ST2, NLAGS, T )
%
% gets                  ST1, ST2        spike times [samples]
%                       NLAGS           maximal time lag [samples]
%                       T               maximal possible spike time [samples]
%
% returns               CCH             histogram of all possible delays,
%                                       corrected for finite durations
%
% calls                 CALC_COR (mex function)
%
% reference             Stark and Abeles, 2009, JNM, pg. 92

% 29-jan-08 ES

% revisions
% 12-dec-19 cleaned up

function cch = calc_cor_flat( st1, st2, NLAGS, T )

if nargin < 4 || isempty( st1 ) || isempty( st2 ) || isempty( NLAGS ) || isempty( T )
    error( 'all arguments are required' )
end

cch1        = calc_cor( st1( st1 <= T - NLAGS ), st2, NLAGS );
cch2        = calc_cor( st2( st2 <= T - NLAGS ), st1, NLAGS );
cch         = [ cch2( 2 * NLAGS + 1 : -1 : NLAGS + 2 ); cch1( NLAGS + 1 : 2 * NLAGS + 1, : ) ];

return

% EOF
