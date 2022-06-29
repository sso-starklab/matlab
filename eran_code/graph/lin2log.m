% lin2log           convert ticks to explicit log scale
%
% CALL              lin2log( ax, base, spacing )
%
% GETS              ax          'x', 'y', or {'xy'}
%                   base        {10}
%                   spacing     {1}
%
% DOES              e.g. if the ticks are -1, 0, 1 and the base is 10, 
%                   the new ticks will be 0.1, 1, 10
% 
% CALLS             nothing

% 24-apr-13 ES

% revisions
% 18-aug-19 cleaned up

function lin2log( ax, base, spacing )

% arguments
nargs                   = nargin;
if nargs < 1 || isempty( ax )
    ax                  = 'xy';
end
if nargs < 2 || isempty( base )
    base                = 10;
end
if nargs < 3 || isempty( spacing )
    spacing             = 1;
end
base                    = base( 1 );

% preps
switch ax
    case 'x'
        lims            = 'xlim';
        fild            = 'xtick';
        lbl             = 'xticklabel';
    case 'y'
        lims            = 'ylim';
        fild            = 'ytick';
        lbl             = 'yticklabel';
    case 'xy'
        lin2log( 'x', base, spacing )
        lin2log( 'y', base, spacing );
        return
end
ticks                   = get( gca, lims );
lo                      = floor( log( base .^ min( ticks ) ) / log( base ) );
hi                      = ceil( log( base .^ max( ticks ) ) / log( base ) );
tick                    = unique( base.^( lo : spacing : hi ) );
logtick                 = log( tick ) / log( base );

% change the labels
set( gca, fild, logtick, lbl, tick )

return

% EOF