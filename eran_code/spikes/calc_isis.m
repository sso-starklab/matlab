% calc_isis                 or cross-isis
%
% call                      y = calc_isis( st1, st2, pad )
%
% gets                      st1         spike train (sparse matrix or list of spike times)
%                           st2         another train
%                           pad         what value to give the first spike
%
% returns                   y           list of ISIs (if st2 not given) or cross-ISIs (if st2 given) 
%
% does                      the ISI is defined as interval between each spike and the preceeding spike
%                           for the first spike, this is undefined, so a NaN (or any value given in 'pad') is assigned
%
%                           if st2 is given, then the cross-ISI is computed
%                           this is defined as the interval between each spike and the preceeding spike from the second train
%                           for any spikes in st1 that occur before the first st2 spike, this is
%                           undefined, so again a NaN is assigned
%
% calls                     nothing
%
% see also                  isisort

% 31-mar-14 ES

% revisions
% 14-oct-19 cleaned up, documented, and checked arrays for being sorted

function y = calc_isis( st1, st2, pad )

y                           = [];
nargs                       = nargin;
if nargs < 2
    st2                     = [];
end
if nargs < 3 || isempty( pad )
    pad                     = NaN;
end
pad                         = pad( 1 );

if issparse( st1 )
    ntrials                 = size( st1, 2 );
    for i                   = 1 : ntrials
        x1                  = find( st1( :, i ) );
        if isempty( x1 )
            continue
        end
        if ~isempty( st2 ) && issparse( st1 ) && isequal( size( st1 ), size( st2 ) )
            x2              = find( st2( :, i ) );
        else
            x2              = [];
        end
        y1                  = calc_isis1( x1, x2, pad );
        y                   = [ y; y1 ];
    end
else
    st1                     = st1( : );
    if isempty( st1 ) 
        return
    end
    if ~isempty( st2 )
        st2                 = st2( : );
    end
    y                       = calc_isis1( st1, st2, pad );
end
                
return

% actual routine:
function y = calc_isis1( st1, st2, pad )
nargs                   = nargin;
if nargs < 2
    st2                 = [];
end
if ~isempty( st1 ) && ~issorted( st1 )
    st1                 = sort( st1 );
end
if isempty( st2 )
    y                   = [ pad; diff( st1 ) ];
    return
end
if ~issorted( st2 )
    st2                 = sort( st2 );
end
i2                      = 1;
y                       = pad * ones( size( st1 ) );
for i1                  = 1 : length( st1 )
    s1                  = st1( i1 );
    if i2 < length( st2 ) && s1 > st2( i2 + 1 )
        i2              = i2 + 1;
    end
    s2                  = st2( i2 );
    if s1 >= s2                                                             % synchrony allowed
        y( i1 )         = s1 - s2;
    end
end

return

% EOF
