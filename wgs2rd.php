<?


function wgs2rd( $lat, $lon ){ // the schrama one
	
	// Fixed constants / coefficients
	$x0      = 155000;
	$y0      = 463000;
	$k       = 0.9999079;
	$bigr    = 6382644.571;
	$m       = 0.003773954;
	$n       = 1.000475857;
	$lambda0 = 0.094032038;
	$phi0    = 0.910296727;
	$l0      = 0.094032038;
	$b0      = 0.909684757;
	$e       = 0.081696831;
	$a       = 6377397.155;
	
	// wgs84 to bessel
	$dphi = $lat - 52;
	$dlam = $lon - 5;
		
	$phicor = ( -96.862 - $dphi * 11.714 - $dlam * 0.125 ) * 0.00001;
	$lamcor = ( $dphi * 0.329 - 37.902 - $dlam * 14.667 ) * 0.00001;

	$phibes = $lat - $phicor;
	$lambes = $lon - $lamcor;

	// bessel to rd
	$phi		= $phibes / 180 * pi();
	$lambda		= $lambes / 180 * pi();
	$qprime		= log( tan( $phi / 2 + pi() / 4 ));  
	$dq			= $e / 2 * log(( $e * sin($phi) + 1 ) / ( 1 - $e * sin( $phi ) ) );
	$q			= $qprime - $dq;
		
	$w			= $n * $q + $m;
	$b			= atan( exp( $w ) ) * 2 - pi() / 2;
	$dl			= $n * ( $lambda - $lambda0 );

	$d_1		= sin( ( $b - $b0 ) / 2 );
	$d_2		= sin( $dl / 2 );
	
	$s2psihalf	= $d_1 * $d_1 + $d_2 * $d_2 * cos( $b ) * cos ( $b0 );
	$cpsihalf	= sqrt( 1 - $s2psihalf );
	$spsihalf	= sqrt( $s2psihalf );
	$tpsihalf	= $spsihalf / $cpsihalf;
	
	$spsi		= $spsihalf * 2 * $cpsihalf;
	$cpsi		= 1 - $s2psihalf * 2;
	
	$sa			= sin( $dl ) * cos( $b ) / $spsi;
	$ca			= ( sin( $b ) - sin( $b0 ) * $cpsi ) / ( cos( $b0 ) * $spsi );
	
	$r			= $k * 2 * $bigr * $tpsihalf;
	$x			= $r * $sa + $x0;
	$y			= $r * $ca + $y0;

	return array('x'=>$x, 'y'=>$y);
}


?>