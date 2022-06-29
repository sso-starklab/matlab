% getWavelet            decomposition of one or more channels
%
% CALL                  [ wave, f, t ] = getWavelet( x, Fs, fMin, fMax, nbins )
%                       [ ..., coh, phases, raw, coi, scale, period, scalef ] = getWaveler( ..., mother, scaling, graphics, nsources)
%
% GETS                  x:      vector - wavelet analysis is applied
%                               multi-column matrix - wavelet analysis is applied to each column separately
%                               two-column matrix - treated as two time series and the csd, coherence is
%                                   computed
%                               3D array, 3rd dimension is 2 - treated as multiple samples from two time series, coherence etc is computed
%
%                       Fs      sampling frequency
%                       fMin, fMax, nbins       parameters for spectral analysis - min/max freqnecy
%                                               (Hz) and number of frequency bins (Morlet is always used so this is
%                                               frequency per se and not just scale)
% 
% OPTIONAL ARGUMENTS
%                       mother          {'morlet'} 
%                       graphics        {0}, flag
%                       scaling         {'var'}
%                       nsources        [], to force no csd 
%
% DOES
% (1) compute the CWT of each segment
% (2) if two signals also compute the cross-spectrum/coherence/phase lag. 
% (3) average over all segments
% (4) plot
%
% RETURNS               wave - the PSD/CSD
%                       f, t - vectors
%                       coh - only for 2-channel input
%                       phases - either for the 1-channel or the phase-difference (from the CSD,
%                               not the smoothed estimate)
% 
% CALLS                 wavelet, smoothwavelet, phaseplot, myjet, colormaps

% 08-oct-12 ES

% revisions
% 15-nov-12 (1) added single-channel phase; cross-spectrum, coherence, and
%               phase difference estimates
%           (2) local plotting of spectrogram + phasogram (single-channel)
%               or coherogram + phase differences (two-channel)
%
% 10-oct-13 (1) mother made a parameter. this changes the input
%               structure... note: check DOG and diff-of-G
% 10-nov-13 (1) scaled output if scaling to Z
% 18-aug-19 cleaned up

% to do: organize input / output handling better, plotting etc
%           also compute mean spectra/coherence/phase 
% also - external coherence computation for multiple segments (i.e. average
%       across trials)

function [ wave, f, t, coh, phases, raw, coi, scale, period, scalef ] = getWavelet( x, Fs, fMin, fMax, nbins, mother, scaling, graphics, nsources )

% arguments
nargs                       = nargin;
if nargs < 5 || isempty( x ) || isempty( Fs ) || isempty( fMin ) || isempty( fMax ) || isempty( nbins )
    error( 'missing arguments' )
end
if isa( x, 'int16' )
    x                       = single( x );
end
if nargin < 6 || isempty( mother )
    mother                  = 'MORLET';
end
if nargin < 7 || isempty( scaling )
    scaling                 = 'var';
end
if nargin < 8 || isempty( graphics )
    if nargout == 0
        graphics            = 1;
    else
        graphics            = 0;
    end
end
if nargin < 9 || isempty( nsources )
    nsources                = [];
end
dt                          = 1 / Fs; 
s0                          = 1 / fMax;
tMax                        = 1 / fMin;
dj                          = log2( tMax / s0 ) / nbins;

% initialize output
wave                        = [];
coh                         = [];
phases                      = [];
raw                         = [];
coi                         = [];
scale                       = [];
period                      = [];
scalef                      = [];

% determine number of signals
er                          = 0;
sx                          = size( x );
switch ndims( x )
    case 2
        switch min( sx )
            case 1
                nsignals        = 1;
                nsegments       = 1;
                nsamples        = max( sx );
                x               = x( : );
            case 2
                if isequal( nsources, 1 )
                    nsignals    = 1;
                    nsegments   = 2;
                else
                    nsignals    = 2;
                    nsegments   = 1;
                    if size( x, 1 ) == 2
                        x       = x';
                    end
                    y           = x( :, 2 );
                    x           = x( :, 1 );
                end
                nsamples        = max( sx );
            otherwise
                nsignals        = 1;
                nsegments       = sx( 2 );
                nsamples        = sx( 1 );
        end
    case 3
        if size( x, 3 ) ~= 2
            er                  = 1;
        else
            nsignals            = 2;
            nsegments           = sx( 2 );
            nsamples            = sx( 1 );
            y                   = x( :, :, 2 );
            x                   = x( :, :, 1 );
        end
    otherwise
        er                      = 1;
end
if er 
    error( 'not supported' )
end

% parameters
if isa( scaling, 'char' ) && strcmp( scaling, 'var' )
    switch nsignals
        case 1
            scalef          = var( x ); % can be a different number for each segment
        case 2
            scalef          = std( x ) .* std( y );
    end
elseif isa( scaling, 'double' ) && size( scaling, 1 ) == ( nbins + 1 )
    scalef                  = scaling; % can be a different number for each frequency
    scaling                 = 'z';
else
    scalef                  = ones( 1, nsegments );
    scaling                 = 'none';
end
flipidx                     = ( nbins + 1 ) : -1 : 1;
t                           = ( 1 : nsamples )' / Fs;

%------------------------------------------------------------------------
% actually compute
%------------------------------------------------------------------------
switch nsignals
    case 1
        xw                  = zeros( nbins + 1, nsamples, nsegments );
        for i               = 1 : nsegments
            [ xw( :, :, i ), period, scale, coi ] = wavelet( x( :, i ), dt, 1, dj, s0, nbins, mother );
        end
        xw                  = xw( flipidx, :, : ); % freq, time, segments
        wave                = abs( xw ) .^ 2;
        phases              = angle( xw );
        raw                 = xw;
    case 2
        xw                  = zeros( nbins + 1, nsamples, nsegments );
        yw                  = xw;
        coh                 = xw;
        for i               = 1 : nsegments
            [ xw( :, :, i ), period, scale, coi ] = wavelet( x( :, i ), dt, 1, dj, s0, nbins, mother );
            [ yw( :, :, i ) ] = wavelet( y( :, i ), dt, 1, dj, s0, nbins, mother );
            % coherence 
            sinv            = 1 ./ ( scale' );
            X               = xw( :, :, i );
            Y               = yw( :, :, i );
            wxy             = X .* conj( Y ); % complex, single trial
            sx              = smoothwavelet( sinv( :, ones( 1, nsamples ) ) .* ( abs( X ) .^ 2 ), dt, period, dj, scale );
            sY              = smoothwavelet( sinv( :, ones( 1, nsamples ) ) .* ( abs( Y ) .^ 2 ), dt, period, dj, scale );
            sWxy            = smoothwavelet( sinv( :, ones( 1, nsamples ) ) .* wxy, dt, period, dj, scale );
            Rsq             = abs( sWxy ) .^ 2 ./ ( sX .* sY );
            coh( :, :, i )  = Rsq( flipidx, : );
        end
        
        % individual channels
        xw                  = xw( flipidx, :, : );
        yw                  = yw( flipidx, :, : );
        raw( :, :, :, 1 )   = xw;
        raw( :, :, :, 2 )   = yw;
        % cross spectrum
        xyw                 = xw .* conj( yw );
        phases              = angle( xyw ); % from the CSD
        wave                = abs( xyw );
        % for multiple trials, one can also compute the trial-averaged coherence and phase lag by:
        %cohTA = mean( abs( xyw ) .^ 2, 3 ) ./ ( mean( abs( xw ) .^ 2, 3 ) .* mean( abs( yw ) .^ 2, 3 ) );     
        %phasesTA = mod( atan2( mean( sin( phases ), 3 ), mean( cos( phases ), 3 ) ), 2 * pi );
end
f                           = 1 ./ period( flipidx );
coi                         = 1 ./ coi;                 % minimum freq to consider at each time point

%------------------------------------------------------------------------
% plot
%------------------------------------------------------------------------
if graphics
    
    % here the scaling is by the signal variance
    
    figure
    
    if nsignals == 1
        nplots              = 1; % PSD
    else
        nplots              = 4; % [ PSD1 CSD; COH PSD2 ]
    end
    
    % scale
    scalename               = '{\sigma}^2';
    pow = zeros( size( wave ) );
    switch scaling
        case { 'var', 'none' }
            scaleres        = 0.25;
            for i           = 1 : nsegments
                pow( :, :, i ) = log2( abs( wave( :, :, i ) / scalef( i ) ) );
            end
        case 'z'
            for i           = 1 : length( f )
                pow( i, :, : ) = ( wave( i, :, : ) - scalef( i, 1 ) ) / scalef( i, 2 );
            end
                scaleres    = ( max( pow( : ) ) - min( pow( : ) ) ) / 100;

    end
    pow                     = mean( pow, 3 );   % average over segments (1 signal)
    levels                  = min( pow(:) ) : scaleres : max( pow( : ) );
    mphases                 = mod( atan2( mean( sin( phases ), 3 ), mean( cos( phases ), 3 ) ), 2 * pi );
    
    % plot
    h1                      = subplot( 1, 1, 1 );
    [ c, h ]                = contourf( t, f, pow, levels );
    set( h, 'linestyle','none')
    ylabel( 'Frequency [Hz]' )
    title( sprintf( '%d segments, %d signals', nsegments, nsignals ) )

    % center color limits around log2(1)=0
    if strcmp( scaling, 'var' )
        clim                = get( gca, 'clim' );
        clim                = [ -1 1 ] * max( clim( 2 ), 3 );
        set( gca, 'clim', clim )
    end
    
    % add the cone of influence
    line( t, coi, 'color', [ 0 0 0 ] );
    hold on
    tt                      = [ t( [ 1 1 ] ) - dt * .5; t; t( [ end end ] ) + dt * .5 ];
    hcoi                    = fill( tt, 1 ./ [ period( [ end 1 ] ) 1./coi period( [ 1 end ] ) ], 'w' );
    set( hcoi, 'alphadatamapping', 'direct', 'facealpha', .5 )
    hold off
    set( h1, 'box', 'off', 'tickdir', 'out' )
    
    % add phase arrows
    if nsignals == 2
        Args.ArrowDensity   = [30 30];
        Args.ArrowSize      = 1;
        Args.ArrowHeadSize  = 1;
        ad                  = mean(Args.ArrowDensity);
        Args.ArrowSize      = Args.ArrowSize*30*.03/ad;
        Args.ArrowHeadSize  = Args.ArrowHeadSize*Args.ArrowSize*220;
        phs_dt              = round(length(t)/Args.ArrowDensity(1));
        tidx                = max(floor(phs_dt/2),1):phs_dt:length(t);
        phs_dp              = round(length(period)/Args.ArrowDensity(2));
        pidx                = fliplr( max(floor(phs_dp/2),1):phs_dp:length(period) );
        phaseplot(t(tidx),f(pidx),2*pi-mphases(pidx,tidx),Args.ArrowSize,Args.ArrowHeadSize);
    end
    
    % add colorbar
    h                       = colorbar;
    subplot( h )
    barylbls                = rats( 2.^ ( get( h, 'ytick' )' ) );
    barylbls( :, all( barylbls == ' ', 1 ) ) = [];
    set( h, 'yticklabel', barylbls );
    title( scalename )
    set( h, 'box', 'off', 'tickdir', 'out' )
    colormap( h1, myjet ) 
    
    figure
    switch nsignals
        case 1
            % plot phases separately
            [ ~, h ]        = contourf( t, f, mphases, 10 );
            set( h, 'linestyle','none')
            xlabel( 'Time [s]' )
            ylabel( 'Frequency [Hz]' )
            set( gca, 'clim', [ 0 2 * pi ] )
            
            % add the cone of influence
            line( t, coi, 'color', [ 0 0 0 ] );
            hold on
            tt = [t([1 1])-dt*.5;t;t([end end])+dt*.5];
            hcoi = fill(tt,1./[period([end 1]) 1./coi period([1 end])],'w');
            set( hcoi, 'alphadatamapping', 'direct', 'facealpha', .5 )
            hold off
            
            h = colorbar;
            subplot( h )
            title( 'Phase (rad)' )
            set( h, 'box', 'off', 'tickdir', 'out' )
            colormap( colormaps( myjet ) )
            
        case 2
            % plot coherence (without phase arrows)
            [ ~, h ]        = contourf( t, f, mean( coh, 3 ), 100 );
            set( h, 'linestyle','none')
            xlabel( 'Time [s]' )
            ylabel( 'Frequency [Hz]' )
            title( sprintf( '%d segments, %d signals', nsegments, nsignals ) )
            set( gca, 'clim', [ 0 1 ] )
            set( h1, 'box', 'off', 'tickdir', 'out' )
            
            % add the cone of influence
            line( t, coi, 'color', [ 0 0 0 ] );
            hold on
            tt              = [t([1 1])-dt*.5;t;t([end end])+dt*.5];
            hcoi            = fill(tt,1./[period([end 1]) 1./coi period([1 end])],'w');
            set( hcoi, 'alphadatamapping', 'direct', 'facealpha', .5 )
            hold off
            
            % add phase plots (only for high coherence values inside the coi)
            aaa             = 2 * pi - mphases;
            aaa( mean( coh, 3 ) < .5 ) = NaN;
            aaa( bsxfun( @lt, f' * ones( 1, nsamples ), coi ) ) = NaN;
            phaseplot( t( tidx ), f( pidx ), aaa( pidx, tidx ), Args.ArrowSize, Args.ArrowHeadSize );
            
            % add colorbar
            h = colorbar;
            subplot( h )
            title( 'Coherence' )
            set( h, 'box', 'off', 'tickdir', 'out' )
            colormap( myjet )
    end
    
    
end

% scaled output
switch scaling
    case 'z'
        for i = 1 : length( f )
            wave( i, :, : ) = ( wave( i, :, : ) - scalef( i, 1 ) ) / scalef( i, 2 );
        end
end

return

% EOF
