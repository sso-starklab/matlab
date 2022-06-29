% stim_make         constructor (empty shell)
%
% s = stim_make
%
% constructor. makes an empty structure with the following fields:
% 
% filebase            string, full path and base
% suffix              string, e.g. 'eeg' or 'dat'
% chan                integer/s; vector if simultaneous structure
% source              string, e.g. 'LED', 'LD', 'DPSS', ..
% voltagerange        [V], e.g. 8
% duration            [sec]
% median              [A2Du] 
% generator           cell array of strings: machine type, date, and a 'concat' flag if generated indirectly by merging multiple structures
% 
% index               n x 1 logical array; 1 if unique, 0 if overlapping with other channels
% types               n x 1 cell array of strings, e.g. 'PULSE', 'SINES', 'SINE',...
% category            n x 3 (or 4) integer matrix: [ type_category value_category duration_category ]
%                           (see stim_categorize for more details)
% durs                n x 1: [sec]
% times               n x 2: [samples]
% slopes              n x 1: [dV/ds]
% vals                n x 1: [V]
% stats               n x 2: mean, SD [V]
% franges             n x 2: freq. at start/end of train, [Hz]; relevant for 'SINES' and 'ZAP'
% plateaus            n x 1; plateau intensity, [V]; relevant for PSINE
%
% see also          stim_check, stim_get

% 18-jan-13 ES

% revisions
% 01-aug-19 exgtended to generate an m by n array of structures

function s = stim_make( m, n )

nargs           = nargin;
if nargs < 1
    m = 1; 
end
if nargs < 2
    n = 1;
end
m               = abs( ceil( m( 1 ) ) );
n               = abs( ceil( n( 1 ) ) );

s.filebase      = '';
s.suffix        = '';
s.chan          = [];
s.source        = '';
s.voltagerange  = [];
s.duration      = [];
s.median        = [];
s.generator     = {};

s.index         = [];            % 1 if unique, {0} if overlapping with other channels
s.types         = [];            % cell array
s.category      = [];         % 2-col integer matrix (according to some partitioning scheme; see stim_categorize)
s.durs          = [];             % [sec]
s.times         = [];            % [samples]
s.slopes        = [];           % [dV/ds]
s.vals          = [];             % [V]
s.stats         = [];            % mean, SD [V]
s.franges       = [];          % chirp/sine only [Hz]
s.plateaus      = [];         % psine only [V]

s               = repmat( s, [ m n ] );

return

% EOF
