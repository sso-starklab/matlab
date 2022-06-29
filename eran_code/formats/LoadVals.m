% LoadVals      load multiple val files in the desired directory
%
% [ Vals, Trigs ] = LoadVals( filebase )
%
% Note: the format of a *.val.* file is
% [ trigger_onset trigger_offset onset_slope offset_slope mean_value ]
% onset/offset are in samples (of the original file, e.g. 20000)
% slopes and mean are in Volts (or could be A2DU)

% 27-jan-12 ES

% 03-sep-12 return a temporally-sorted array
% 17-aug-19 cleaned up

function [ Vals, Trigs ] = LoadVals( filebase )

valbase             = sprintf( '%s.val*', filebase );
valdir              = fileparts( filebase );
valfnames           = dir( valbase );

Vals                = [];
Trigs               = [];
% accumulate
for i               = 1 : length( valfnames )
    valfname        = sprintf( '%s/%s', valdir, valfnames( i ).name );
    val             = LoadVal( valfname );
    if valfnames( i ).name( end - 2 ) == 't'
        chan        = str2num( valfnames( i ).name( end - 1 : end ) );
    else
        chan        = str2num( valfnames( i ).name( end - 2 : end ) );
    end
    Vals            = [ Vals; val ];
    Trigs           = [ Trigs; chan * ones( size( val, 1 ), 1 ) ];
end
% sort
if ~isempty( Vals )
    [ Vals, sidx ]  = sortrows( Vals, 1 );
    Trigs           = Trigs( sidx );
end

return

% EOF
