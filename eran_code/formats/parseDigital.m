% parseDigital          parse a 16-bit channel into 16 one-bit events
%
% call                  mat = parseDigital( filebase )
%
% optional arguments, given as name/value pairs:
%
%                       chan            otherwise will extract from par
%                       nchans          otherwise will extract from par
%                       spkFs           otherwise will extract from par
%
%                       shiftbits       {2^15}
%                       precision       {'int16'}   If 'uint16', set shiftbits to 0
%                       nbits           {16}
%                       Overwrite       {-2}
% 
% returns               mat, 3 column matrix: [ time bitnum diff ]
%                           time:       samples 
%                           bitnum:     1..16
%                           diff:       low-to-high, 1
%                                       high-to-low, -1
%
% calls                 if chan/nchans/spkFs not given:     LoadXml, get_stimchans
%                       otherwise:                          nothing
%
% Notes:
% 1. This routine is slow, can take almost 1 minute for 1 hour recording
% 2. At the time of writing, bit 16 was constantly on, so although data
%           were recorded as int16, routine was called with shiftbits=0
%   Once corrected, shiftbits should maintain its default value of 2^15
%   Using the bit extractor without the shift would be correct only for uint16
 
% 01-sep-18 ES

% revisions
% 02-sep-18 added saving of *evt.d## files
% 30-oct-18 (1) modified typo in call to ParseArgPairs (was nbits twice, 
%               first one should be spkFs)
%           (2) added verbal report of number of blocks
% 24-dec-20 added option for digstr for over 100 channels

function mat = parseDigital( filebase, varargin )

%----------------------------------------------------------------------%
% preps
%----------------------------------------------------------------------%

% constants
blocksize       = 1e6;      % [samples/channel]

% arguments
nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ chan, nchans, spkFs ...
    , shiftbits, precision, nbits, Overwrite ...
    ] = ParseArgPairs(...
    { 'chan', 'nchans', 'spkFs' ...
    , 'shiftbits', 'precision', 'nbits', 'Overwrite' ...
    }...
    , { [], [],  [] ...
    , 2^15, 'int16', 16, -2 ...
    }...
    , varargin{ : } );
if isempty( chan ) || isempty( nchans ) || isempty( spkFs )
    par             = LoadXml( filebase );
end
if isempty( chan )
    chan        = get_stimchans( par, [], 'digitalin' );
end
if isempty( nchans )
    nchans      = par.nChannels;
end
if isempty( spkFs )
    spkFs       = par.SampleRate;
end
if str2double( precision( end - 1 : end ) ) ~= nbits 
    error( 'input mismatch' )
end

% parameters
eval( sprintf( 'v = %s( 1 );', precision ) );
vi              = whos( 'v' ); 
nbytes          = vi.bytes;
fname           = [ filebase '.dat' ];
bitnum          = 1 : nbits;

% load if requested and exists
if chan < 10
    digstr      = sprintf( '00%d', chan );
elseif chan < 100
    digstr      = sprintf( '0%d', chan );
elseif chan < 1000
    digstr      = sprintf( '%d', chan );
end
digfname        = sprintf( '%s.dig.%s', filebase, digstr );
dexists         = exist( digfname, 'file' );
if Overwrite < 0 && dexists
    fprintf( 1, '%s: Loading existing parsed digital data for %s\n', upper( mfilename ), filebase )
    load( digfname, '-mat', 'mat' );
    return
else
    fprintf( 1, '%s: Parsing digital data for %s, channel #%d\n', upper( mfilename ), filebase, chan )
end       

%----------------------------------------------------------------------%
% parse
%----------------------------------------------------------------------%

% partition into blocks
info            = dir( fname );
nsamples        = info.bytes / nbytes / nchans;
nblocks         = ceil( nsamples / blocksize );
blocks          = [ 1 : blocksize : blocksize * nblocks; blocksize : blocksize : blocksize * nblocks ]';
blocks( nblocks, 2 ) = nsamples;
fprintf( '\t\t%d blocks\t', nblocks )

% one-sample overlap to prevent missing transitions at block transision
blocks( 2 : nblocks, 1 ) = blocks( 2 : nblocks, 1 ) - 1;

% go over blocks and parse
mat             = [];
for i           = 1 : nblocks
    
    % extract the bits
    boff        = ( blocks( i, 1 ) - 1 ) * nbytes * nchans;
    bsize       = ( diff( blocks( i, : ) ) + 1 );
    m           = memmapfile( fname, 'Format', precision, 'Offset', boff, 'Repeat', bsize * nchans );
    d           = m.Data( chan : nchans : bsize * nchans );
    clear m
    dmat        = rem( floor( single( d + shiftbits ) * 2.^( 1 - bitnum ) ), 2 );

    % detect the transitions
    bmat        = [];
    for j = 1 : nbits
        dvec    = diff( dmat( :, j ) );
        hi      = blocks( i, 1 ) + find( dvec == 1 );
        lo      = blocks( i, 1 ) + find( dvec == -1 );
        hiblock = [ hi ones( length( hi ), 1 ) * [ j  1 ] ];
        loblock = [ lo ones( length( lo ), 1 ) * [ j -1 ] ];
        bmat    = [ bmat; hiblock; loblock ];
    end    
    bmat        = sortrows( bmat, 1 );

    % accumulate
    mat         = [ mat; bmat ];
    
    % report
    if ~mod( i, 10 )
        fprintf( 1, '.' )
    end
    
end
mat             = sortrows( mat, 1 );
fprintf( 1, 'Done.\n' )

%----------------------------------------------------------------------%
% save
%----------------------------------------------------------------------%

if Overwrite == 1 || ( Overwrite ~= -1 && ~dexists )
    
    % save matlab format binary for all bits
    fprintf( 1, '%s: Saving parsed digital data for %s\n', upper( mfilename ), filebase )
    save( digfname, 'mat' )
    
    % save a text format event file for each bit (Neuroscope-compatible) 
    fprintf( 1, 'saving event files for ' )
    dignums             = unique( mat( :, 2 ) );
    for i               = 1 : length( dignums )
        dignum          = dignums( i );
        fprintf( 'bit#%d ', dignum )
        if dignum < 10
            digstr      = sprintf( '0%d', dignum );
        elseif i < 100
            digstr      = sprintf( '%d', dignum );
        end
        evtfname        = sprintf( '%s.evt.d%s', filebase, digstr );
        
        didx            = mat( :, 2 ) == dignum;
        evtime          = mat( didx, 1 ) / spkFs * 1000;                        % convert to ms
        evlabel         = mat( didx, 3 );                                       % assign labels: 1: rising edge; -1: falling edge
        fid             = fopen( evtfname, 'wt' );
        fprintf( fid, '%10.2f  %2.0f\n', [ evtime( : )'; evlabel( : )' ] );     % full resolution <=100h of recording...
        fclose( fid );
    end
    fprintf( 1, 'done!\n' )
    
end

return

% EOF
