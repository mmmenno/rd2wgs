<?


function rd2wgs ($x, $y){
	
    // Calculate WGS84 coördinates
    $dX = ($x - 155000) * pow(10, - 5);
    $dY = ($y - 463000) * pow(10, - 5);
    $SomN = (3235.65389 * $dY) + (- 32.58297 * pow($dX, 2)) + (- 0.2475 *
         pow($dY, 2)) + (- 0.84978 * pow($dX, 2) *
         $dY) + (- 0.0655 * pow($dY, 3)) + (- 0.01709 *
         pow($dX, 2) * pow($dY, 2)) + (- 0.00738 *
         $dX) + (0.0053 * pow($dX, 4)) + (- 0.00039 *
         pow($dX, 2) * pow($dY, 3)) + (0.00033 * pow(
            $dX, 4) * $dY) + (- 0.00012 *
         $dX * $dY);
    $SomE = (5260.52916 * $dX) + (105.94684 * $dX * $dY) + (2.45656 *
         $dX * pow($dY, 2)) + (- 0.81885 * pow(
            $dX, 3)) + (0.05594 *
         $dX * pow($dY, 3)) + (- 0.05607 * pow(
            $dX, 3) * $dY) + (0.01199 *
         $dY) + (- 0.00256 * pow($dX, 3) * pow(
            $dY, 2)) + (0.00128 *
         $dX * pow($dY, 4)) + (0.00022 * pow($dY,
            2)) + (- 0.00022 * pow(
            $dX, 2)) + (0.00026 *
         pow($dX, 5));
 
    $Latitude = 52.15517 + ($SomN / 3600);
    $Longitude = 5.387206 + ($SomE / 3600);
 
    return array(
        'latitude' => $Latitude ,
        'longitude' => $Longitude
		);
}



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