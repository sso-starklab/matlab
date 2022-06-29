% INLINE_STATS      inlines for basic statistical tests.
%
% call              [ BINP, PP, CHI2P, CHI2XY ] = INLINE_STATS
%
% returns           inline functions; compute the probability of getting
%
%                   BINP        obs or more when exp are expected (cumulative binomial; exact, one-sided)
%                   PP          obs or more when exp are expected (cumulative poisson; one sided)
%                   CHI2P       obs when p * tot are expected (chi2 test; approximate; two-sided)
%                   CHI2XY      x out of totx and y out of toty the same by chance (chi2 test)
%
% calls             nothing; the inlines call BINOMIAL, POISSCDF, and CHI2_TEST
%
% called by         PLOT_VPM, SUA_SUM (SUA_PACK).

% 29-may-04 ES

% revisions
% 17-jun-04 binomial -> binocdf
% 17-jul-04 CHI2P corrected

function [ binp, pp, chi2p, chi2xy ] = inline_stats;

%binp = inline( '1 - binomial( obs - 1, tot, p )', 'obs', 'tot', 'p' );
binp = inline( '1 - binocdf( obs - 1, tot, p )', 'obs', 'tot', 'p' );
pp = inline( '1 - poisscdf( obs - 1, exp )', 'obs', 'exp' );
%chi2p = inline( 'chi2_test( [ obs tot; ceil( tot * p ) tot] )', 'obs', 'tot', 'p' );
% chi2p = inline( 'chi2_test( [ obs tot - obs; ceil( tot * p ) tot - ceil( tot * p ) ] )', 'obs', 'tot', 'p' );
% chi2xy = inline( 'chi2_test( [ x totx - x; y toty - y ] )', 'x', 'totx', 'y', 'toty' );
chi2p = inline( 'gtest( [ obs tot - obs; ceil( tot * p ) tot - ceil( tot * p ) ] )', 'obs', 'tot', 'p' );
chi2xy = inline( 'gtest( [ x totx - x; y toty - y ] )', 'x', 'totx', 'y', 'toty' );

return