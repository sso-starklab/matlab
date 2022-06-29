% a single trial by digital bits of DTDT
% bits 0-15

% 1) 8 OR 9       - PSs box

% 8 is on R arm? 

% 2) 2            - spotter after PWs
% 
%    trial left:
%       3) 11     - PS left arm
%       4) 6      - sol left
%
%   right trial:
%       3) 10     - PS right arm
%       4) 7      - sol right
%
% 5) 12-15        - status bits
%
% 6) 9 OR 8       - PSs box

function dtdt_trial

if ismac 
    filebase = '/Volumes/slab1/mA335/dat/mA335_49/mA335_49';
elseif ispc || islinux
    filename                        = 'mA335_49'; filebase = filebase_lookup(filename,2);
end

load( [ filebase '.mov.mat' ] );

dchan = get_stimchans( filebase, [], 'digitalin' );
load( [ filebase '.dig.' num3str( dchan ) ], '-mat' );
[ aa, bb ] = uhist( mat( :, 2 ) );
% [ bb' aa' ]
% 
%            2         586
%            3         270
%            5        4029
%            6      205475
%            7         100
%            8          98
%            9         152
%           10         174
%           11         230
%           12         340

% simple events:
% 3     CA
% 7     L_Sol
% 8     R_Sol
% 9     R_box 
% 10    L_box
% 11    R_Arm
% 12    L_Arm


%
% sequence:
% 11 := choice -> (11) -> 9 := new trial -> (9,10) -> 
% 12 := choice -> (12) -> 10 := new trial -> (9,10) ->

% 11 + 8 := correct
% 10 + 7 := correct

% L_Sol                           = mat( ismember( mat( :, 2 : 3 ), [ 7 1 ], 'rows' ), 1 );
% R_Sol                           = mat( ismember( mat( :, 2 : 3 ), [ 8 1 ], 'rows' ), 1 );
% R_Box                           = mat( ismember( mat( :, 2 : 3 ), [ 9 1 ], 'rows' ), 1 );
% L_Box                           = mat( ismember( mat( :, 2 : 3 ), [ 10 1 ], 'rows' ), 1 );
% R_Arm                           = mat( ismember( mat( :, 2 : 3 ), [ 11 1 ], 'rows' ), 1 );
% L_Arm                           = mat( ismember( mat( :, 2 : 3 ), [ 12 1 ], 'rows' ), 1 );
% 
% m1                              = [ L_Sol ones( length( L_Sol ), 1 ) * 7 ];
% m2                              = [ R_Sol ones( length( R_Sol ), 1 ) * 8 ];
% m3                              = [ R_Box ones( length( R_Box ), 1 ) * 9 ];
% m4                              = [ L_Box ones( length( L_Box ), 1 ) * 10 ];
% m5                              = [ R_Arm ones( length( R_Arm ), 1 ) * 11 ];
% m6                              = [ L_Arm ones( length( L_Arm ), 1 ) * 12 ];
% 
% m                               = [ m1; m2; m3; m4; m5; m6 ];
% m                               = sortrows( m, 1 );

%----------------------------------------------------------------
% Simple event definitions
CA      = 3;
L_Sol   = 7;
R_Sol   = 8;
R_box   = 9;
L_box   = 10;
R_Arm   = 11;
L_Arm   = 12;

%val                             = 1;                                        % we are looking for onset (low->high), which are 1
val                             = -1;                                       % but this failed in this dataset; thus, we used -1
m                               = mat( mat( :, 3 ) == val & ismember( mat( :, 2 ) , [ 3 7 : 12 ] ), 1 : 2 );
%m                               = mat( ismember( mat( :, 2 ) , [ 3 7 : 12 ] ), : );


%----------------------------------------------------------------
% Composite events definition
% Turn_R            :=  first 11 after 9 or 10 (time of 11)
% Turn_L            :=  first 12 after 9 or 10 (time of 12)

% Rtrn_R            :=  first 9 after Turn_R
% Rtrn_L            :=  first 10 after Turn_L

% Trl_srt           := Rtrn_R || Rtrn_L
% Correct_R         := Turn_R && 8
% Correct_L         := Turn_L && 7

%----------------------------------------------------------------
% Turn_R:
pmat                            = m( ismember( m( :, 2 ), [ R_box L_box ] ), : );
tmat                            = m( m( :, 2 ) == R_Arm, : );
% pmat = pmat( 1 : 20, : )
% tmat = tmat( 1 : 20, : )
nt                              = size( tmat, 1 );
np                              = size( pmat, 1 );
Turn_R                          = zeros( min(nt,np), 2 );
Rtrn_R                          = zeros( min(nt,np), 2 );
i                               = 1;
j                               = 1;
k                               = 0;
p                               = pmat( i, 1 );
t                               = tmat( j, 1 );
while i <= np && j <= nt
    while p < t
        i                       = i + 1;
        if i > np
            break
        end
        p                       = pmat( i, 1 );
    end
    k                           = k + 1;
    Turn_R( k, : )              = tmat( j, : );
    if i<=np
        Rtrn_R (k, : )              = pmat( i, : );
    end
    while t < p
        j                       = j + 1;
        if j > nt
            break
        end
        t                       = tmat( j, 1 );
    end
end
Turn_R( k + 1 : end, : )         = [];
Rtrn_R( k + 1 : end, : )         = [];

%----------------------------------------------------------------
% Turn_L:
pmat                            = m( ismember( m( :, 2 ), [ R_box L_box ] ), : );
tmat                            = m( m( :, 2 ) == L_Arm, : );
nt                              = size( tmat, 1 );
np                              = size( pmat, 1 );
Turn_L                          = zeros( min(nt,np), 2 );
Rtrn_L                          = zeros( min(nt,np), 2 );
i                               = 1;
j                               = 1;
k                               = 0;
p                               = pmat( i, 1 );
t                               = tmat( j, 1 );
while i <= np && j <= nt
    while p < t
        i                       = i + 1;
        if i > np
            break
        end
        p                       = pmat( i, 1 );
    end
    k                           = k + 1;
    Turn_L( k, : )              = tmat( j, : );
    if i<=np
        Rtrn_L (k, : )              = pmat( i, : );
    end
    while t < p
        j                       = j + 1;
        if j > nt
            break
        end
        t                       = tmat( j, 1 );
    end
end
Turn_L( k + 1 : end, : )         = [];
Rtrn_L( k + 1 : end, : )         = [];

%----------------------------------------------------------------

%----------------------------------------------------------------
% Composite events definition
% Turn_R            :=  first R_Arm after R_box or L_box (time of R_Arm)
% Turn_L            :=  first L_Arm after R_box or L_box (time of L_Arm)

% Rtrn_R            :=  first R_box after Turn_R
% Rtrn_L            :=  first L_box after Turn_L




% tmat                = parse( find( m( :, 2 ) == L_Arm ) );
% Turn_L              = m( tmat( :, 1 ), : );
% 
% tmat                = parse( find( m( :, 2 ) == R_Arm ) );
% Turn_R              = m( tmat( :, 1 ), : );

% tmat                = m( m( :, 2 ) == R_box, : );
% for i               = 1 : length( Turn_R )
%     j               = find( tmat( :, 1 ) > Turn_R( i ), 1, 'first' );
%     Rtrn_R( i, : )  = tmat( j, : );
% end

return

% EOF

