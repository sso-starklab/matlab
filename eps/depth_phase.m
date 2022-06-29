colors_NPB              = [ 75 54 145; 173 20 87; 42 101 176] / 255; % Nunit, Punit, Biphasic
colors_IP_light         = [ 96 173 94; 156 77 204 ] / 255; % INT, PYR

didx                = ~isnan( tD.depth );
    
tD1                = struct_select( tD, didx );

isbip           = ~isnan(tD1.bpi) & ~isinf( tD1.bpi );
ispos           = tD1.extremum > 0 & ~isbip;
ispyr           = tD1.shankclu (:,3) == 1;

% Slabels         = struct_select( tD1, ~ispos & ~isbip);
% plabels         = p(~ispos & ~isbip);
% ridx            = isnan(plabels);
% dlabels         = d(~ispos & ~isbip);
% plabels (ridx)  = [];
% dlabels (ridx)  = [];
% Slabels         = struct_select( Slabels, ~isnan(plabels));
% labels          = Slabels.shankclu (:,3);
%     f                   = tD1.geo_fwhm;
%     w                   = 1 ./ tD1.fmax * 1000; % width
    t                   = tD1.t2p;               % t2p
    amp                 = abs( tD1.extremum );
    r                   = tD1.FR;            % firing rate
    a                   = tD1.ach_com;          % burstiness
    p                   = tD1.prefPhs (:,1);
    fnames              = tD1.session_tag;
    
    % compute depth in um based on probe resolution
    ufnames             = unique( tD1.session_tag );
    nufnames            = length( ufnames );
    tmpcell             = strfind( ufnames, 'mC41' ); % presently supporting only one animal
    % to support multiple animals (strings), must expand the code
    tmpvec              = NaN( nufnames, 1 );
    for i               = 1 : nufnames
        if isempty( tmpcell{ i } )
            tmpvec( i ) = 0;
        else
            tmpvec( i ) = tmpcell{ i };
        end
    end
    tmpvec              = logical( tmpvec );
    hridx               = ismember( fnames, ufnames( tmpvec ) );
    d( hridx )          = tD1.depth( hridx ) * 15;         % [um]
    d( ~hridx )         = tD1.depth( ~hridx ) * 20;        % [um]
    
    % determine spatial bins:
    binsize             = 10; 
    minVal              = floor( min( d ) / binsize ) * binsize;
    maxVal              = ceil( max( d ) / binsize ) * binsize;
    rside               = ( 0 + binsize / 2 ) : binsize : ( maxVal + binsize / 2 );
    lside               = ( 0 + binsize / 2 ) : binsize : ( -minVal + binsize / 2 );
    edges               = [ fliplr( -lside ) rside ];
    binC                = ( edges( 1 : end - 1 )  + edges( 2 : end ) )' / 2;
    
    figure
    hold on, 
%     ph( 1 ) = plot( p( ~ispos & ~isbip & ispyr ), d( ~ispos & ~isbip & ispyr ), '.' ); 
%     ph( 2 ) = plot( p( ~ispos & ~isbip & ~ispyr ), d( ~ispos & ~isbip & ~ispyr ), '.' ); 
    ph( 3 ) = plot( p( ispos ), d( ispos ), '.' ); 
    set( ph, 'MarkerSize', 6 )
    ph( 4 ) = plot( p( isbip ), d( isbip ), '.' ); 
    set( ph, 'MarkerSize', 10 )
    set( ph( 1 ), 'color', colors_IP_light( 2, : ) )    
    set( ph( 2 ), 'color', colors_IP_light( 1, : ) )    
    set( ph( 3 ), 'color', colors_NPB( 2, : ) )
    set( ph( 4 ), 'color', colors_NPB( 3, : ) )
    xlabel( 'Phase lock' )
    ylabel( 'Depth [\mum]' )
    set( gca, 'tickdir', 'out', 'box', 'off' )
    alines( 0, 'y', 'color', [ 0 0 0 ], 'linestyle', '--' );   
    set( gca, 'xtick', 0 : pi / 2 : 2 * pi )
    set( gca, 'XTickLabel', round( ( 0 : pi / 2 : 2 * pi ) / 0.1 ) * 0.1 )