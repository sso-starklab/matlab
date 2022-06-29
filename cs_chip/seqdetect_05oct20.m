% seqdetect             detect sequence of events according to a rule
%
% call                  [ seqs, smat ] = seqdetect( mats, seqmat )
%
% gets                  mats        cell array of matrices of ranges.
%                                       one matrix per sequence element
%
%                       seqmat      matrix describing the rule. 
%                                       each row corresponds to one sequence
%                                       element. 
%                                   syntax of each row: [ mindur maxdur minlag maxlag ]
%                                   lags are from onset of element i to onset of element i+1
%
% returns               seqs        [ nseqs by 2 ] matrix of all ranges
%                                   each row is a range of one sequence
%
%                       smat        [ nseqs x nelements x 4 ] matrix
%                                   each row contains: 
%                                       -the range of one event
%                                       -event number
%                                       -the interval from the previous event
%
% example:              seqmat = [ 5 15 25 30; 10 20 NaN NaN ];
% 
%                       this corresponds to a 5-15 sample event in mats{ 1 }, followed by a 10-20
%                       sample event in mats{ 2 }. The mats{ 2 } event must lag the mats{ 1 }
%                       event by 25-30 samples.
%
% calls                 sortranges              (sets)
%                       tempmatch               (ssp)

% 04-oct-20 ES

function [ seqs, smat ] = seqdetect( mats, seqmat, dt )


%------------------------------------------------------------------
% check arguments
seqs                            = [];
smat                            = [];

%------------------------------------------------------------------
% check arguments
nargs                           = nargin;
if nargs < 2 || isempty( mats ) || isempty( seqmat )
    error( 'missing arguments' )
end
if nargs < 3 || isempty( dt )
    dt                          = 0;
end
if ~isa( mats, 'cell' )
    error( 'input type mismatch: mats must be cell array' )
end
mats                            = mats( : );
if ~isa( seqmat, 'numeric' ) || ~ismatrix( seqmat )
    error( 'seqmat must be a matric of numbers' )
end
ne                              = length( mats );
if size( seqmat, 1 ) ~= ne
    error( 'input size mismatch: seqmat must have same number of rows as the number of matrices in mats' )
end
for i                           = 1 : ne
    mat                         = mats{ i };
    if ~isa( mat, 'numeric' ) || ~ismatrix( mat ) || size( mat, 2 ) ~= 2 
        fprintf( 1, 'element %d of mats must be a matrix of ranges', i )
        return
    end
    if ~isequal( sortranges( mat ), mat )
        fprintf( 1, 'element %d of mats must be a matrix of non-overlapping ranges', i )
        return
    end
end

%------------------------------------------------------------------
% (1) prune elements of each matrix by durations
nev                             = zeros( ne, 1 );
for i                           = 1 : ne
    mat                         = mats{ i };
    durs                        = diff( mat, [], 2 ) + dt;
    drange                      = seqmat( i, 1 : 2 );
    didx                        = durs >= drange( 1 ) & durs <= drange( 2 );
    mats{ i }                   = mat( didx, : );
    nev( i )                    = sum( didx );
end

%------------------------------------------------------------------
% (2) combine all events
ei                              = cumsum( nev );
tev                             = ei( ne );
si                              = [ 1; ei( 1 : ne - 1 ) + 1 ];
emat                            = zeros( tev, 4 );
for i                           = 1 : ne
    ridx                        = si( i ) : ei( i );
    emat( ridx, 1 : 3 )         = [ mats{ i } ones( nev( i ), 1 ) * i ];
end
emat                            = sortrows( emat, 1 );
difs                            = [ 0; emat( 2 : tev, 1 ) - emat( 1 : tev - 1, 2 ) ];
emat( :, 4 )                    = difs;

%------------------------------------------------------------------
% (3) detect sequences, regardless of lags
tmpl                            = ( 1 : ne )';
c                               = tempmatch( tmpl, emat( :, 3 ) );
c                               = abs( c );
tmpj                            = tmpl( [ 1 : ne - 1 ne - 1 ] );
z1                              = ( tmpl - mean( tmpl ) ) / std( tmpl );
z2                              = ( tmpj - mean( tmpj ) ) / std( tmpj );
th                              = sum( z1 .* z2 ) / ( ne - 1 );
th                              = ( th + 1 ) / 2;
tidx                            = find( c >= th );
ns                              = length( tidx );
m1                              = tidx * ones( 1, ne );
m2                              = ones( ns, 1 ) * ( 0 : ne - 1 );
m                               = ( m1 + m2 )';
smat0                           = emat( m( : ), : );

%------------------------------------------------------------------
% (4) keep the sequences with the proper lags
dlags                           = [ 0 inf; seqmat( 1 : ne - 1, 3 : 4 ) ];
lmat                            = repmat( dlags, [ ns 1 ] );
tf                              = smat0( :, 4 ) >= lmat( :, 1 ) & smat0( :, 4 ) <= lmat( :, 2 );
tf                              = reshape( tf, [ ns ne ] );
tf                              = sum( tf, 2 ) == ne;
sidx                            = 1 + ( find( tf ) - 1 ) * ne;
nseqs                           = length( sidx );
m1                              = sidx * ones( 1, ne );
m2                              = ones( nseqs, 1 ) * ( 0 : ne - 1 );
m                               = ( m1 + m2 )';
smat                            = smat0( m( : ), : );
seqs                            = [ smat0( sidx, 1 ) smat0( sidx + ne - 1, 2 ) ];

return

% EOF

%------------------------------------------------------------------
% 04-oct-20 example run: 

% load the data (durations and onset times)
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs42.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs43.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs45.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/durs47.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims42.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims43.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims45.mat')
load('/Users/eranstark/Documents/da/cs_chip/2020_10_04/tims47.mat')

% convert durations to samples, and create range matrices
mat42                           = tims42 + [ zeros( length( durs42 ), 1 ) round( durs42 * 20000 ) - 1 ];
mat43                           = tims43 + [ zeros( length( durs43 ), 1 ) round( durs43 * 20000 ) - 1 ];
mat45                           = tims45 + [ zeros( length( durs45 ), 1 ) round( durs45 * 20000 ) - 1 ];
mat47                           = tims47 + [ zeros( length( durs47 ), 1 ) round( durs47 * 20000 ) - 1 ];

% combine into cell array at set the rule
mats                            = { mat42, mat45, mat43, mat47 };
seqmat                          = [ 7 12 29 31; 7 12 29 31; 7 12 29 31; 7 12 NaN NaN ] * 20;

% call the routine
[ seqs, smat ]                  = seqdetect( mats, seqmat, 1 );
