<?php
/**
 * rdnaptrans.php
 * 
 * This program is officially validated and may carry the following trademark: RDNAPTRANS™2018
 * It is a combination of two php functionthat performs RD NAP -> ETRS89 transformations or v.v.
 * Coordinate transformation to and from Stelsel van de Rijksdriehoeksmeting and Normaal Amsterdams Peil.
 * - etrs2rd(): Transformation from ETRS89 to RD.
 * - rd2etrs():Transformation from RD to ETRS89 (called the inverse transformation).
 * 
 * This implementation follows the flow as described in the RDNAPTRANS2018-document found at www.nsgi.nl
 * 
 * As a kind of bonus the transformation from ETRS89 to ITRS2014 and v.v. is also included.
 * 
 * I would like to thank Jochem Lesparre from the Netherlands Kadaster/NSGI for his kind help with finding the bugs.  
 */

$GLOBALS['correctionCache']=array();


/**
 * prefillCache
 * 
 * Reads the complete correctiongrid and the height-matrix into the global-variable.
 *
 * @return void
 */
function prefillCache() {
    if (count($GLOBALS['correctionCache']) == 144782) { return; }

    //ETRS89_lat_(deg)	ETRS89_lon_(deg)	NAP_quasi_geoid_height_above_ETRS89_ellipsoid_(m)
    //50.0000	2.0000	44.6078

    foreach(file("nlgeo2018.txt") as $line) {
        if (substr($line,0,1) != 'E') { //skips the first line
            $arr=explode("\t",$line);

            $lat=(float)$arr[0];      
            $lon=(float)$arr[1];      
            $h=(float)$arr[2];
            
            $index=($lat*10000).($lon*10000);
            $GLOBALS['correctionCache'][$index]=array();
            $GLOBALS['correctionCache'][$index]['Height']=$h;
        }
    }

    //RD_lat_(deg)	RD_lon_(deg)	lat_corr_(deg)	lon_corr_(deg)
    //50.0000	 2.0000	+0.000000000	+0.000000000

    foreach(file("rdcorr2018.txt") as $line) {
        if (substr($line,0,1) != 'E') { //skips the first line
            $arr=explode("\t",$line);

            $lat=(float)$arr[0];      
            $lon=(float)$arr[1];      
            $latcor=(float)$arr[2];
            $loncor=(float)$arr[3];

            $index=($lat*10000).($lon*10000);
            $GLOBALS['correctionCache'][$index]['LatCor']=$latcor;
            $GLOBALS['correctionCache'][$index]['LonCor']=$loncor;
        }
    }
}
  
/**
 * etrs2rd
 * 
 * Performs the calculation from ETRD89 to RD.
 *
 * @param  array $etrs   Array with the following keys:
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 *      h   height in meters (optional)
 * @return array $rd    Array with the following keys:
 *      x   x-value of the location in meters
 *      y   y-value of the location in meters
 *      h   height in meters (no height given returns NaN)
 */
function etrs2rd($etrs) {
    prefillCache();

    $grs80=array('a'=>6378137,'f'=>1/298.257222101,'h'=>43);
    $bessel1841=array('a'=>6377397.155,'f'=>1/299.1528128);
    $etrs2rd=array('tx'=>-565.7346,'ty'=>-50.4058,'tz'=>-465.2895,'rx'=>-0.395023,'ry'=>0.330776,'rz'=>-1.876073,'s'=>-4.07242E-6);
    $rdmap=array('latc'=>deg_to_rad1(52.156160556),'lonc'=>deg_to_rad1(5.387638889),'k'=>0.9999079,'x0'=>155000,'y0'=>463000, 'ellips'=>$bessel1841);

    $radian=deg_to_rad($etrs);

    $etrs89=ellips_to_geoc($radian, $grs80);

    $rdxyz=helmert($etrs89, $etrs2rd);

    $rd_pseudo_bessel_rad=geoc_to_ellips($rdxyz, $bessel1841);

    $corrected_deg=RD_correction($rd_pseudo_bessel_rad);
    
    $corrected_rad=deg_to_rad($corrected_deg);

    $rd=ellips_to_rd($corrected_rad, $rdmap);

    if ($etrs['h']!==null) {
        $rd['z']=NAP_correction($etrs);
    } else {
        $rd['z']=NAN;
    }
 
    return $rd;
}

/**
 * rd2etrs
 * 
 * Performs the calculation from RD to ETRD89.
 *
 * @param array $rd    Array with the following keys:
 *      x   x-value of the location in meters
 *      y   y-value of the location in meters
 *      h   height in meters (optional)
 * @return  array $etrs   Array with the following keys:
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 *      h   height in meters (no height given returns NaN)
 */
function rd2etrs($rd) {
    prefillCache();
    $grs80=array('a'=>6378137,'f'=>1/298.257222101,'h'=>43);
    $bessel1841=array('a'=>6377397.155,'f'=>1/299.1528128,'h'=>0);
    $rdmap=array('latc'=>deg_to_rad1(52.156160556),'lonc'=>deg_to_rad1(5.387638889),'k'=>0.9999079,'x0'=>155000,'y0'=>463000, 'ellips'=>$bessel1841);
    $rd2etrs=array('tx'=>565.7381,'ty'=>50.4018,'tz'=>465.2904,'rx'=>0.395026,'ry'=>-0.330772,'rz'=>1.876074,'s'=>4.07242E-6);

    $rd_real_bessel_rad=rd_to_ellips($rd, $rdmap);

    $rd_real_bessel_deg=rad_to_deg($rd_real_bessel_rad);

    $rd_pseudo_bessel_deg=inverse_RD_Correction($rd_real_bessel_deg);

    $rd_pseudo_bessel_rad=deg_to_rad($rd_pseudo_bessel_deg);

    $rdxyz=ellips_to_geoc($rd_pseudo_bessel_rad, $bessel1841);
    
    $etrsxyz=helmert($rdxyz, $rd2etrs);
    
    $etrs_rad=geoc_to_ellips($etrsxyz, $grs80);

    $etrs_deg=rad_to_deg($etrs_rad);

    $etrs_deg['h']=$rd['h'];
    $etrs_deg['h']=NAP_Correction($etrs_deg,false);
    
    return $etrs_deg;
}

/**
 * ITRS2ETRS89
 *
 * Performs the calculation from ITRS to ETRS89 (or the inverse).
 *
 * @param  array $pos   Array with the following keys:
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 *      h   height in meters (optional, if not given h=0)
 * 
 * @param  boolean $inverse  If true inverses the calculation
 * @param  mixed $date Date to be used in the calcution. 
 *      Could be a DateTime-object or a string containg a date 
 * @return array
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 *      h   height in meters
 */
function ITRS2ETRS89($pos, $inverse=false, $date=null) {
    $grs80=array('a'=>6378137,'f'=>1/298.257222101,'h'=>43);

    $itrf2014=array('tx'=>0.0547,'ty'=>0.0522,'tz'=>-0.0741,'rx'=>-0.001701,'ry'=>-0.010290,'rz'=>+0.016632,'s'=>+0.00000000212);
    $itrf2014_rate=array('tx'=>0.00010,'ty'=>0.00010,'tz'=>-0.00190,'rx'=>-0.0000810,'ry'=>-0.0004900,'rz'=>+0.0007920,'s'=>+0.000000000110);

    $pos=deg_to_rad($pos);

    $phi=$pos['lat'];
    $labda=$pos['lon'];
    $h=0;
    if (isset($pos['h'])) {
        $h=$pos['h'];
        $grs80['h']=$h;
    }

    $phi1 = $phi;
    
    $geoc=ellips_to_geoc($pos, $grs80);

    if (gettype($date)=='string') {
        $date=date_create($date);
    }

    if ($date === null || !($date instanceof DateTime)) {
        $date=new datetime();
    }

    $t_epoch=2010;
    $t=$date->format("Y")+round($date->format('z')/366,2);

    $wgs2itrf=helmert_apply_rate($itrf2014, $itrf2014_rate, $t-$t_epoch);

    if ($inverse) {
        $wgs2itrf=inverse_helmert_param($wgs2itrf);
    }

    $helmert=helmert($geoc, $wgs2itrf);

    $grs80=array('a'=>6378137,'f'=>1/298.257222101,'h'=>43);
    $ellips=geoc_to_ellips($helmert, $grs80);

    return rad_to_deg($ellips);
} 

/**
 * inverse_RD_Correction
 *
 * @param  array Geographic ellipsoidal real Bessel in degrees 
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 * @return array Geographic ellipsoidal pseudo Bessel in decimal degrees.
 *      lat latitude in decimal degrees
 *      lon longitude in decimal degrees
 */
function inverse_RD_Correction($pos) {
    $phi = $pos['lat'];        // in degrees;
    $labda =  $pos['lon'];

    $phi_min = 50; 
    $phi_max = 56;
    $labda_min = 2;
    $labda_max = 8;
    $phi_delta = 0.0125;
    $labda_delta = 0.02;

    $phinorm = ($phi-$phi_min)/$phi_delta;
    $labdanorm = ($labda-$labda_min)/$labda_delta;
    $nlabda = 1 + (($labda_max-$labda_min)/$labda_delta);

    if ($phi >= $phi_min && $phi <= $phi_max && $labda >= $labda_min && $labda <= $labda_max ) {
  
        $NorthWest=getRdCorrection($phi, $labda, 1, 0);
        $nw_phi=$NorthWest['LatCor'];
        $nw_labda=$NorthWest['LonCor'];

        $NorthEast=getRdCorrection($phi, $labda, 1, 1);
        $ne_phi=$NorthEast['LatCor'];
        $ne_labda=$NorthEast['LonCor'];

        $SouthWest=getRdCorrection($phi, $labda, 0, 0);
        $sw_phi=$SouthWest['LatCor'];
        $sw_labda=$SouthWest['LonCor'];

        $SouthEast=getRdCorrection($phi, $labda, 0, 1);
        $se_phi=$SouthEast['LatCor'];
        $se_labda=$SouthEast['LonCor'];

        $RDcorrLat = ($phinorm - floor($phinorm)) * (($nw_phi  *(floor($labdanorm) + 1 - $labdanorm)) + $ne_phi  *($labdanorm - floor($labdanorm))) + (floor($phinorm) + 1 - $phinorm) * (($sw_phi  *(floor($labdanorm) + 1 - $labdanorm)) + $se_phi  *($labdanorm - floor($labdanorm)));
        $RDcorrLon = ($phinorm - floor($phinorm)) * (($nw_labda*(floor($labdanorm) + 1 - $labdanorm)) + $ne_labda*($labdanorm - floor($labdanorm))) + (floor($phinorm) + 1 - $phinorm) * (($sw_labda*(floor($labdanorm) + 1 - $labdanorm)) + $se_labda*($labdanorm - floor($labdanorm)));

        $xperc=($labdanorm - floor($labdanorm));
        $ncor=$xperc*$ne_phi+(1-$xperc)*$nw_phi;
        $scor=$xperc*$se_phi+(1-$xperc)*$sw_phi;
        $yperc=($phinorm - floor($phinorm));
        $RDcorrLat=$yperc*$ncor+(1-$yperc)*$scor;

        $xperc=($labdanorm - floor($labdanorm));
        $ncor=$xperc*$ne_labda+(1-$xperc)*$nw_labda;
        $scor=$xperc*$se_labda+(1-$xperc)*$sw_labda;
        $yperc=($phinorm - floor($phinorm));
        $RDcorrLon=$yperc*$ncor+(1-$yperc)*$scor;
    } else {
        $RDcorrLat = 0;
        $RDcorrLon = 0;
    }
    $phi1 = $phi + $RDcorrLat;
    $labda1 = $labda + $RDcorrLon;

    return array('lat'=>$phi1,'lon'=>$labda1);
} // inverse_RD_Correction;
  
/**
 * rd_to_ellips
 *
 * @param  array RD X coordinate in meter
 * @param  array Geographic ellipsoidal real Bessel latitude in radian
 *      lat latitude in radian
 *      lon longitude in radian
 * @return void
 */
function rd_to_ellips($rd,$map) {
    // Inverse oblique stereographic conformal projection from the RD projection plane to a sphere;
    $const_90_degrees = pi()/2;              // 90 degrees in radian;
    $const_180_degrees = pi();               // 180 degrees in radian;
    $const_360_degrees = 2*pi();             // 360 degrees in radian;
    $x=$rd['x'];
    $y=$rd['y'];


    // Get the parameters of RD Map projection and Bessel 1841 ellipsoid parameter;
    $phi0 = $map['latc']; //symget('gmv_phi0_amersfoort'); 
    $labda0 = $map['lonc']; //symget('gmv_labda0_amersfoort'); 
    $k = $map['k']; //symget('gmv_k_amersfoort');
    $x0 = $map['x0']; //symget('gmv_x0_amersfoort');
    $y0 = $map['y0']; //symget('gmv_y0_amersfoort');
    $a = $map['ellips']['a']; //symget('gmv_B1841_a');
    $f = $map['ellips']['f']; //symget('gmv_B1841_f');
    $epsilon = 0.00000000002; //symget('gmv_epsilon_B1841_threshold');

    // Start with derived parameter calculation of the RD map projection;
    $e =  sqrt($f*(2-$f));
    $RN = $a / sqrt(1 - ($e**2*(sin($phi0)**2)));
    $RM = ($RN*(1 - $e**2))/(1-($e**2*(sin($phi0)**2)));
    $R_sphere = sqrt($RM*$RN);
    $PHI0_C = atan((sqrt($RM)/sqrt($RN))*tan($phi0)); 
    $LABDA0_C = $labda0;

    // Inverse oblique stereographic projection of coordinates on plane to coordinates on sphere;
    $r_distance = sqrt(($x - $x0)**2 + ($y - $y0)**2);
    $sin_alpha = ($x - $x0)/$r_distance;
    $cos_alpha = ($y - $y0)/$r_distance;
    $psi =  2 * atan($r_distance/(2*$k*$R_sphere));
    if ($x != $x0 || $y != $y0 ) {
        $Xnorm = cos($PHI0_C)*cos($psi) - $cos_alpha*sin($PHI0_C)*sin($psi);
        $Ynorm = $sin_alpha*sin($psi);
        $Znorm = $cos_alpha*cos($PHI0_C)*sin($psi) + sin($PHI0_C)*cos($psi);
    } else {
        $Xnorm = cos($PHI0_C);
        $Ynorm = 0;
        $Znorm = sin($PHI0_C);
    }
    $PHI_C = asin($Znorm);
    if ($Xnorm > 0 ) { $LABDA_C = $LABDA0_C + atan($Ynorm/$Xnorm);
    } elseif ($Xnorm < 0 && $x >= x0 ) { $LABDA_C = $LABDA0_C + atan($Ynorm/$Xnorm) + $const_180_degrees;
    } elseif ($Xnorm < 0 && $x < x0 ) { $LABDA_C = $LABDA0_C + atan($Ynorm/$Xnorm) - $const_180_degrees;
    } elseif ($Xnorm == 0 && $x > x0 ) { $LABDA_C = $LABDA0_C + $const_90_degrees;
    } elseif ($Xnorm == 0 && $x < x0 ) { $LABDA_C = $LABDA0_C - $const_90_degrees;
    } elseif ($Xnorm == 0 && $x == x0 ) { $LABDA_C = $LABDA0_C;
    } 

    /*
    Projection from sphere to ellipsoid: The second step of the inverse RD map projection is an inverse
    Gauss conformal projection from the sphere to the Bessel ellipsoid to obtain Bessel coordinates.
    Start with remaining derived parameter calculation of the RD map projection.
    */
    $q0 = log(tan(($phi0 + $const_90_degrees)/2)) - ($e/2) * (log((1 + $e*sin($phi0))/(1 - $e*sin($phi0))));
    $w0 = log(tan(($PHI0_C + $const_90_degrees)/2));
    $n = sqrt(1 + (($e**2*(cos($phi0)**4))/(1 - $e**2)));
    $m = $w0 - $n*$q0;

    // Inverse Gauss conformal projection of coordinates on sphere to coordinates on ellipsoid;
    $w = log(tan(($PHI_C + $const_90_degrees)/2));
    $q = ($w - $m )/$n;    
 
    $phi1 = $PHI_C;
    $i = 0;
    do {
        $phi = $phi1;
        if ($PHI_C > -1*$const_90_degrees && $PHI_C < $const_90_degrees ) { 
            $phi1 = 2*atan(exp($q + ($e/2)*log((1 + $e*sin($phi))/(1 - $e*sin($phi))))) - $const_90_degrees;
        } else {$phi1 = $PHI_C;}
        $i++;
    } while (abs($phi1-$phi) >= $epsilon);
    
    // The latitute has been calculated, now calculate the longitude;
    $labda_n = (($LABDA_C - $LABDA0_C)/$n) + $labda0;
    $labda = $labda_n + $const_360_degrees*floor(($const_180_degrees - $labda_n)/$const_360_degrees);

    return array('lat'=>$phi1, 'lon'=>$labda);
} // rd_to_ellips

/**
 * getRdCorrection
 *
 * @param  mixed $lat
 * @param  mixed $lon
 * @param  mixed $offsetx
 * @param  mixed $offsety
 * @return void
 */
function getRdCorrection($lat, $lon, $offsetx=null, $offsety=null) {
    $LatInt=floor($lat/0.0125)*125;
    $LonInt=floor($lon/0.02)*200;

    if ($offsetx !== null && $offsetx !== 0) {$LatInt+=125;}
    if ($offsety !== null && $offsety !== 0) {$LonInt+=200;}

    if ($LatInt <500000 || $LatInt>560000 ) return array('LatCor'=>0,'LatCor'=>0, 'Height'=>null);
    if ($LonInt < 20000 || $LonInt> 80000 ) return array('LatCor'=>0,'LatCor'=>0, 'Height'=>null);


    if ( array_key_exists($LatInt.$LonInt, $GLOBALS['correctionCache'])) {
        return $GLOBALS['correctionCache'][$LatInt.$LonInt];
    }
}

/**
 * RD_correction
 *
 * @param  array $pos    Pseudo Bessel RD in radian.
 * @return array The geographic ellipsoidal real Bessel latitude in decimal degrees
 */
function RD_correction($pos) {
    $posdeg=rad_to_deg($pos);

    $poscompare=$posdeg;

    $phi_min = 50; 
    $phi_max = 56; 
    $labda_min = 2; 
    $labda_max = 8; 
    $phi_delta = 0.0125; 
    $labda_delta = 0.02; 
    $epsilon = 1.0E-9;
    $phi_threshold= false;
    $labda_threshold = false;
    $currentpos1=$posdeg;

    do {

        $phi = $currentpos1['lat'];
        $labda = $currentpos1['lon'];

        $phinorm = ($phi-$phi_min)/$phi_delta;
        $labdanorm = ($labda-$labda_min)/$labda_delta;
        $nlabda = 1 + (($labda_max-$labda_min)/$labda_delta);
        
        if ($phi > $phi_min && $phi <= $phi_max && $labda >= $labda_min && $labda <= $labda_max ) {
            $inside_bound_correction_grid = 1;

            $NorthWest=getRdCorrection($phi, $labda, 1, 0);
            $nw_phi=$NorthWest['LatCor'];
            $nw_labda=$NorthWest['LonCor'];

            $NorthEast=getRdCorrection($phi, $labda, 1, 1);
            $ne_phi=$NorthEast['LatCor'];
            $ne_labda=$NorthEast['LonCor'];

            $SouthWest=getRdCorrection($phi, $labda, 0, 0);
            $sw_phi=$SouthWest['LatCor'];
            $sw_labda=$SouthWest['LonCor'];

            $SouthEast=getRdCorrection($phi, $labda, 0, 1);
            $se_phi=$SouthEast['LatCor'];
            $se_labda=$SouthEast['LonCor'];
        } else {
            $inside_bound_correction_grid = 0;
        }

        if (! $phi_threshold) {
            if ($inside_bound_correction_grid == 1) {
                $RDcorrLat = ($phinorm - floor($phinorm)) * (($nw_phi*(floor($labdanorm) + 1 - $labdanorm)) + $ne_phi*($labdanorm - floor($labdanorm))) + (floor($phinorm) + 1 - $phinorm) * (($sw_phi*(floor($labdanorm) + 1 - $labdanorm)) + $se_phi*($labdanorm - floor($labdanorm)));
            } else {
                $RDcorrLat = 0;
            }

            $currentpos1['lat'] = $poscompare['lat'] - $RDcorrLat;

            if (abs($currentpos1['lat']-$phi) < $epsilon) {
                $phi_threshold = true;
            }
        }

        if (! $labda_threshold) {
            if ($inside_bound_correction_grid == 1 ) { 
                $RDcorrLon = ($phinorm - floor($phinorm)) * (($nw_labda*(floor($labdanorm) + 1 - $labdanorm)) + $ne_labda*($labdanorm - floor($labdanorm))) + (floor($phinorm) + 1 - $phinorm) * (($sw_labda*(floor($labdanorm) + 1 - $labdanorm)) + $se_labda*($labdanorm - floor($labdanorm)));
            } else {
                $RDcorrLon = 0;
            }

            $currentpos1['lon'] = $poscompare['lon'] - $RDcorrLon;

            if (abs($currentpos1['lon']-$labda) < $epsilon) {
                $labda_threshold = true;
            }
        }

    } while ($phi_threshold === false || $labda_threshold === false);

    return array('lat'=>$currentpos1['lat'],'lon'=>$currentpos1['lon']);
} // RD_correction


/**
 * ellips_to_rd
 *
 * @param  array $pos
 * @param  array $map
 * @return array
 */
function ellips_to_rd($coor, $map) {
    $e = sqrt($map['ellips']['f']*(2-$map['ellips']['f']));
    $q0 = log(tan(($map['latc'] + pi()/2)/2)) - ($e/2) * (log((1 + $e*sin($map['latc']))/(1 - $e*sin($map['latc']))));
    $Rn = $map['ellips']['a'] / sqrt(1 - ($e*$e*(sin($map['latc'])**2)));
    $RM = ($Rn*(1 - $e*$e))/(1-($e*$e*(sin($map['latc'])**2)));
    $R_sphere = sqrt($RM*$Rn);
    $PHI0_C = atan((sqrt($RM)/sqrt($Rn))*tan($map['latc'])); 
    $LABDA0_C = $map['lonc'];
    $w0 = log(tan(($PHI0_C + pi()/2)/2));
    $n = sqrt(1 + (($e*$e*pow(cos($map['latc']),4))/(1 - $e*$e)));
    $m = $w0 - $n*$q0;

    $q = log(tan(($coor['lat'] + pi()/2)/2)) - ($e/2) * (log((1 + $e*sin($coor['lat']))/(1 - $e*sin($coor['lat']))));
    $w = $n*$q + $m;
    $PHI_C =  2*atan(exp($w)) - pi()/2;
    $LABDA_C = $LABDA0_C + $n*($coor['lon'] - $map['lonc']);

    $sin_psi_2 = sqrt(pow(sin(($PHI_C - $PHI0_C)/2),2) + (pow(sin(($LABDA_C - $LABDA0_C)/2),2)*cos($PHI_C)*cos($PHI0_C)));
    $cos_psi_2 = sqrt(1 - pow($sin_psi_2,2));
    $tan_psi_2 = $sin_psi_2/$cos_psi_2;
    $sin_alpha = (sin($LABDA_C - $LABDA0_C)*cos($PHI_C))/(2*$sin_psi_2*$cos_psi_2);
    $cos_alpha = (sin($PHI_C) - sin($PHI0_C) + 2*sin($PHI0_C)*($sin_psi_2*$sin_psi_2))/(2*cos($PHI0_C)*$sin_psi_2*$cos_psi_2);
    $r_distance = 2*$map['k']*$R_sphere*$tan_psi_2;

    if ($PHI_C == $PHI0_C && $LABDA_C == $LABDA0_C) {;
        $RD_x = $map['x0'];
        $RD_y = $map['y0'];
    } else {
        if (($PHI_C != $PHI0_C || $LABDA_C != $LABDA0_C) && ($PHI_C != -1*$PHI0_C || $LABDA_C != pi()-$LABDA0_C)) {
            $RD_x = $r_distance*$sin_alpha + $map['x0'];
            $RD_y = $r_distance*$cos_alpha + $map['y0'];
        } else {
            $RD_x = 0;
            $RD_y = 0;
        }
    }

    return array('x'=>$RD_x, 'y'=>$RD_y, 'z'=>43);
}

/**
 * NAP_Correction
 *
 * @param  array $pos
 * @param  bool $ETRS89toRD
 * @return array
 */
function NAP_Correction($pos, $ETRS89toRD=true) {

    $phi=$pos['lat'];
    $labda=$pos['lon'];
    if (isset($pos['h'])) {
        $height=$pos['h'];
    } else {
        return NAN;
    }

    $phi_min = 50; 
    $phi_max = 56; 
    $labda_min = 2; 
    $labda_max = 8; 
    $phi_delta = 0.0125; 
    $labda_delta = 0.02; 

    $phinorm = ($phi-$phi_min)/$phi_delta;
    $labdanorm = ($labda-$labda_min)/$labda_delta;
    $nlabda = 1 + (($labda_max-$labda_min)/$labda_delta);
  
    // Rounding needed for the test validation;
    if (round($phi*10000000)/10000000 >= $phi_min && round($phi*10000000)/10000000 <= $phi_max && round($labda*10000000)/10000000 >= $labda_min && round($labda*10000000)/10000000 <= $labda_max ) {
  
        $NorthWest=getRdCorrection($phi, $labda, 1, 0);
        $nw_height=$NorthWest['Height'];

        $NorthEast=getRdCorrection($phi, $labda, 1, 1);
        $ne_height=$NorthEast['Height'];

        $SouthWest=getRdCorrection($phi, $labda, 0, 0);
        $sw_height=$SouthWest['Height'];

        $SouthEast=getRdCorrection($phi, $labda, 0, 1);
        $se_height=$SouthEast['Height'];

        $etrs89_quasi_height = ($phinorm - floor($phinorm)) * (($nw_height*(floor($labdanorm) + 1 - $labdanorm)) + $ne_height*($labdanorm - floor($labdanorm))) + (floor($phinorm) + 1 - $phinorm) * (($sw_height*(floor($labdanorm) + 1 - $labdanorm)) + $se_height*($labdanorm - floor($labdanorm)));
    } else {
        return NAN;
    }
    if ($ETRS89toRD) {
        return $height - $etrs89_quasi_height;
    } else {
        return $height + $etrs89_quasi_height;
    }
} //NAP_Correction;

/**
 * deg_to_rad1
 *
 * @param  float $deg  Location in degrees
 * @return float Location in radian
 */
function deg_to_rad1($deg) {
    return ($deg * pi())/180;
} //deg_to_rad1

/**
 * deg_to_rad
 *
 * @param  array Location in degrees
 * @return array location in radian
 */
function deg_to_rad($pos) {
    $phi=deg_to_rad1($pos['lat']);
    $labda=deg_to_rad1($pos['lon']);

    if (isset($pos['h'])) {
        return array('lat'=>$phi,'lon'=>$labda,'h'=>$pos['h']);
    } else {
        return array('lat'=>$phi,'lon'=>$labda);
    }
} // deg_to_rad

/**
 * rad_to_deg
 *
 * @param  array $pos  Location in radian
 * @return array Location in degrees
 */
function rad_to_deg($pos) {
    
    $phi=($pos['lat']*180)/pi();
    $labda=($pos['lon']*180)/pi();

    if (isset($pos['h'])) {
        return array('lat'=>$phi,'lon'=>$labda,'h'=>$pos['h']);
    } else {
        return array('lat'=>$phi,'lon'=>$labda);
    }
} // rad_to_deg

/**
 * ellips_to_geoc
 *
 * @param  array $pos    Coordinates in elipsiodal coordinates
 * @param  mixed $param  Parameters of the ellipsoid
 * @return array Coordinates in geocantric position
 */
function ellips_to_geoc($pos, $param) {
    $phi=$pos['lat'];
    $labda=$pos['lon'];

    $a=$param['a'];
    $f=$param['f'];
    $h=$param['h'];

    $e2=$f*(2-$f);

    $Rn=$a/sqrt(1-$e2*pow(sin($phi),2));

    $X=($Rn + $h)*cos($phi)*cos($labda);
    $Y=($Rn + $h)*cos($phi)*sin($labda);
    $Z=($Rn*(1-$e2) + $h)*sin($phi);

    return array('x'=>$X, 'y'=>$Y, 'z'=>$Z);
} // ellips_to_geoc

/**
 * geoc_to_ellips
 *
 * @param  array $pos    Location in geocantric position 
 * @param  array $param  Parameters of the ellipsoid
 * @return array elipsiodal coordinates
 */
function geoc_to_ellips($pos, $param) {
    $e2=$param['f']*(2-$param['f']);
    $epsilon=2E-10;
    $phi1=45;

    do {
        $phi=$phi1;

        $Rn=$param['a']/sqrt(1-$e2*pow(sin($phi1),2));

        if ($pos['x'] > 0) {
            $phi1 = atan(($pos['z'] + $e2*$Rn*sin($phi))/sqrt(pow($pos['x'],2) + pow($pos['y'],2)));}
        else {
            if (round($pos['x'],11) == 0 && round($pos['y'],11) == 0 && round($pos['z'],11) >= 0 ) {
                $phi1 = pi()/2;
            } else {
                $phi1 = -1*pi()/2;
            }
        }                                                                                                // -90 degrees;
    } while (!(abs($phi1-$phi) < $epsilon));

    if ($pos['x'] > 0) { $labda = atan($pos['y']/$pos['x']);
    } elseif ($pos['x'] < 0 and $pos['y'] >= 0 ) { $labda = atan($pos['y']/$pos['x']) + pi();// +180 degrees;
    } elseif ($pos['x'] < 0 and $pos['y'] < 0 ) { $labda = atan($pos['y']/$pos['x']) - pi();// -180 degrees;
    } elseif ($pos['x'] == 0 and $pos['y'] > 0 ) { $labda = pi()/2;          //  +90 degrees;
    } elseif ($pos['x'] == 0 and $pos['y'] < 0 ) { $labda = -1*pi()/2;       //  -90 degrees;
    } elseif ($pos['x'] == 0 and $pos['y'] == 0 ) { $labda = 0;}

    $h=sqrt(pow($pos['x'],2)+pow($pos['y'],2))*cos($phi1)+$pos['z']*sin($phi1)-$param['a']*sqrt(1-$e2*pow(sin($phi1),2));

    return array('lat'=>$phi1, 'lon'=>$labda, 'h'=>$h);
} // geoc_to_ellips

/**
 * helmert_apply_rate
 *
 * @param  array $source  Transformation parameters
 * @param  array $rate    Change rate of the transformation parameters
 * @param  float $t
 * @return array Transformation parameters with the change applies
 */
function helmert_apply_rate($source, $rate, $t) {
    $result=array();

    $result['tx']=$source['tx']+$rate['tx']*$t;
    $result['ty']=$source['ty']+$rate['ty']*$t;
    $result['tz']=$source['tz']+$rate['tz']*$t;
    $result['rx']=$source['rx']+$rate['rx']*$t;
    $result['ry']=$source['ry']+$rate['ry']*$t;
    $result['rz']=$source['rz']+$rate['rz']*$t;
    $result['s']=$source['s']+$rate['s']*$t;
    
    return $result;
}

/**
 * inverse_helmert_param
 *
 * @param  array $param  Transformation parameters
 * @return array Transformation parameters, for the inverse function
 */
function inverse_helmert_param($param) {
    $result=array();
    $result['tx']=-$param['tx'];
    $result['ty']=-$param['ty'];
    $result['tz']=-$param['tz'];
    $result['rx']=-$param['rx'];
    $result['ry']=-$param['ry'];
    $result['rz']=-$param['rz'];
    $result['s'] =-$param['s'];
    
    return $result;
}

/**
 * helmert
 * 
 * Performs a 3D Helmert transformation on geocentric coordinates. 
 *
 * @param  array $source  Geocentric position
 * @param  array $param   Transformation parameters
 * @return array Transformed geocentric position
 */
function helmert($source, $param) {
    $param['rx']=deg_to_rad1($param['rx']/3600);
    $param['ry']=deg_to_rad1($param['ry']/3600);
    $param['rz']=deg_to_rad1($param['rz']/3600);

    $s=1+$param['s'];
    $R11=cos($param['rz'])*cos($param['ry']);
    $R12=cos($param['rz'])*sin($param['ry'])*sin($param['rx'])+sin($param['rz'])*cos($param['rx']);
    $R13=-cos($param['rz'])*sin($param['ry'])*cos($param['rx'])+sin($param['rz'])*sin($param['rx']);

    $R21=-sin($param['rz'])*cos($param['ry']);
    $R22=-sin($param['rz'])*sin($param['ry'])*sin($param['rx'])+cos($param['rz'])*cos($param['rx']);
    $R23=sin($param['rz'])*sin($param['ry'])*cos($param['rx'])+cos($param['rz'])*sin($param['rx']);

    $R31=sin($param['ry']);
    $R32=-cos($param['ry'])*sin($param['rx']);
    $R33=cos($param['ry'])*cos($param['rx']);


    $X=$s*($R11*$source['x']+$R12*$source['y']+$R13*$source['z'])+$param['tx'];
    $Y=$s*($R21*$source['x']+$R22*$source['y']+$R23*$source['z'])+$param['ty'];
    $Z=$s*($R31*$source['x']+$R32*$source['y']+$R33*$source['z'])+$param['tz'];

    return array('x'=>$X, 'y'=>$Y, 'z'=>$Z);
} // helmert
?>