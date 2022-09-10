<?php

if (isset($_REQUEST['m'])) {
    include 'rdnaptrans.php';
    if ($_REQUEST['m']==1) {
        $rd=array('x'=>108360.8790,'y'=>415757.2745,'z'=>258.0057);
        printf ('RD: X: %10.4f, Y: %10.4f, Z: %10.4f<br/>', $rd['x'], $rd['y'], $rd['z']);
        $etrs89=rd2etrs($rd);
        printf ('ETRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $etrs89['lat'], $etrs89['lon'], $etrs89['h']);

    } elseif ($_REQUEST['m']==2) {
        $etrs89=array('lat'=>53.473095072,'lon'=>6.886681267,'h'=>290.4741);
        printf ('ETRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $etrs89['lat'], $etrs89['lon'], $etrs89['h']);
        $rd=etrs2rd($etrs89);
        printf ('RD: X: %10.4f, Y: %10.4f, Z: %10.4f<br/>', $rd['x'], $rd['y'], $rd['z']);
    } elseif ($_REQUEST['m']==3) {
        SelfValidation();
    } elseif ($_REQUEST['m']==4) {
        ValidateETRS2RD();
    } elseif ($_REQUEST['m']==5) {
        ValidateRD2ETRS();
    } elseif ($_REQUEST['m']==6) {
        $etrs89=array('lat'=>51.999994814, 'lon'=>4.999992333, 'h'=>42.9821);
        printf ('ETRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $etrs89['lat'], $etrs89['lon'], $etrs89['h']);
        $itrs=ITRS2ETRS89($etrs89, true, '2020-10-1');
        printf ('ITRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $itrs['lat'], $itrs['lon'], $itrs['h']);
    } elseif ($_REQUEST['m']==7) {
        $itrs=array('lat'=>52,'lon'=>5,'h'=>43);
        printf ('ITRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $itrs['lat'], $itrs['lon'], $itrs['h']);
        $etrs89=ITRS2ETRS89($itrs, false, '2020-10-1');
        printf ('ETRS: lat: %13.9f, lon: %13.9f, H: %10.4f<br/>', $etrs89['lat'], $etrs89['lon'], $etrs89['h']);
    }
} else {
?>

<a href="?m=1">Convert 1 demo-coordinate from RD to ETRS89</a></br>
<a href="?m=2">Convert 1 demo-coordinate from ETRS89 to RD</a></br>
<a href="?m=3">Perform selfvalidation</a></br>
<a href="?m=4">Create validationfile for ETRS2RD</a></br>
<a href="?m=5">Create validationfile for RD2ETRS</a></br>
<a href="?m=6">Convert 1 demo-coordinate from ETRS89 to ITRS</a></br>
<a href="?m=7">Convert 1 demo-coordinate from ITRS to ETRS89</a></br>
<?php
}

function ValidateRD2ETRS() {
    header('Content-Description: File Transfer');
    header('Content-Disposition: attachment; filename=etrstext.txt');
    header('Expires: 0');
    header('Cache-Control: must-revalidate');
    header('Pragma: public');
    header("Content-Type: text/plain");
    foreach(file("002_RDNAP.txt") as $line) {
        //point_id  x_coordinate  y_coordinate height
        //10020000   127767.3607   839078.2344  38.7032
        //01234567890123456789012345678901234567890123456789
        $name=substr($line,0,8);
        $x=(double)substr($line,10,12);
        $y=(double)substr($line,24,12);
        $h=(double)substr($line,37,8);

        //point_id  latitude        longitude    height
        //20020000  52.754767993    3.006916444  302.1601
        //01234567890123456789012345678901234567890123456789
        if (substr($name,0,1)=='p') {
            echo 'point_id  latitude        longitude    height'."\r\n";
        } else {
            $etrs=rd2etrs(array('x'=>$x,'y'=>$y,'h'=>$h));
            if (is_nan($etrs['h'])) {
                printf("%8s  %12.9f   %12.9f       NaN"."\r\n",$name,$etrs['lat'],$etrs['lon']);
            } else {
                printf("%8s  %12.9f   %12.9f %9.4f"."\r\n",$name,$etrs['lat'],$etrs['lon'],$etrs['h']);
            }
        }
    }
    return;
}


function ValidateETRS2RD() {
    header('Content-Description: File Transfer');
    header('Content-Disposition: attachment; filename=rdtext.txt');
    header('Expires: 0');
    header('Cache-Control: must-revalidate');
    header('Pragma: public');
    header("Content-Type: text/plain");
    foreach(file("002_ETRS89.txt") as $line) {
        $name=substr($line,0,8);
        $lat=(double)substr($line,9,13);
        $lon=(double)substr($line,25,13);
        $h=(double)substr($line,39,9);

        if (substr($name,0,1)=='p') {
            echo 'point_id  x_coordinate  y_coordinate height'."\r\n";
        } else {
            $rd=etrs2rd(array('lat'=>$lat,'lon'=>$lon,'h'=>$h));
            if (is_nan($rd['z'])) {
                printf("%8s  %12.4f  %12.4f        NaN"."\r\n",$name,$rd['x'],$rd['y']);
            } else {
                printf("%8s  %12.4f  %12.4f  %9.4f"."\r\n",$name,$rd['x'],$rd['y'],$rd['z']);
            }
        }

    }
    return;
}

/**
 * SelfValidation
 *
 * @return void
 */
function SelfValidation() {
    $i=1;
    $fout=0;
    echo '<table border="1px black solid">';
    foreach(file("Z001_ETRS89andRDNAP.txt") as $line) {
        $arr=explode("\t",$line);

        $rd=etrs2rd(array('lat'=>$arr[1],'lon'=>$arr[2],'h'=>(double)$arr[3]));
        $etrs=rd2etrs(array('x'=>$arr[4],'y'=>$arr[5],'h'=>(double)$arr[6]));
        $chk=false;
        $chk=$chk | (abs($arr[1]-$etrs['lat'])>0.00000001);
        $chk=$chk | (abs($arr[2]-$etrs['lon'])>0.00000001);
        $chk=$chk | (abs((double)$arr[3]-$etrs['h'])>0.001 && substr($arr[6],0,1)!='*');
        $chk=$chk | (abs($arr[4]-$rd['x'])>0.0001);
        $chk=$chk | (abs($arr[5]-$rd['y'])>0.0001);
        $chk=$chk | (abs((double)$arr[6]-$rd['z'])>0.0001);
        
        if ($chk) {
            $fout++;
            echo '<tr>';
            printf ('<td>%s X</td>', $i);

            printf ('<td>%s</td>', $arr[0]);
            if (abs($arr[1]-$etrs['lat'])<0.00000001) {
                printf ('<td>%12.9f</td>', $etrs['lat']);
            } elseif (abs($arr[1]-$etrs['lat'])<0.0000001) {
                printf ('<td><font color="green">%12.9f</font></br>%12.9f</td>', $etrs['lat'],$arr[1]);
            } else {
                printf ('<td><font color="red">%12.9f</font></br>%12.9f</td>', $etrs['lat'],$arr[1]);
            }
            if (abs($arr[2]-$etrs['lon'])<0.00000001) {
                printf ('<td>%12.9f</td>', $etrs['lon']);
            } elseif (abs($arr[2]-$etrs['lon'])<0.0000001) {
                printf ('<td><font color="green">%12.9f</font></br>%12.9f</td>', $etrs['lon'],$arr[2]);
            } else {
                printf ('<td><font color="red">%12.9f</font></br>%12.9f</td>', $etrs['lon'],$arr[2]);
            }

            if (abs((double)$arr[3]-$etrs['h'])<0.0001 || substr($arr[6],0,1)=='*') {
                printf ('<td>%8.4f</td>', $etrs['h']);
            } elseif (abs((double)$arr[3]-$etrs['h'])<0.001) {
                printf ('<td><font color="green">%8.4f</font></br>%8.4f</td>', $etrs['h'],$arr[3]);
            } else {
                printf ('<td><font color="red">%8.4f</font></br>%8.4f</td>', $etrs['h'],$arr[3]);
            }

            if (abs($arr[4]-$rd['x'])<0.0001) {
                printf ('<td>%8.4f</td>', $rd['x']);
            } elseif (abs($arr[4]-$rd['x'])<0.001) {
                printf ('<td><font color="green">%8.4f</font></br>%8.4f</td>', $rd['x'],$arr[4]);
            } else {
                printf ('<td><font color="red">%8.4f</font></br>%8.4f</td>', $rd['x'],$arr[4]);
            }

            if (abs($arr[5]-$rd['y'])<0.0001) {
                printf ('<td>%8.4f</td>', $rd['y']);
            } elseif (abs($arr[5]-$rd['y'])<0.001) {
                printf ('<td><font color="green">%8.4f</font></br>%8.4f</td>', $rd['y'],$arr[5]);
            } else {
                printf ('<td><font color="red">%8.4f</font></br>%8.4f</td>', $rd['y'],$arr[5]);
            }

            if (abs((double)$arr[6]-$rd['z'])<0.0001 || substr($arr[6],0,1)=='*') {
                printf ('<td>%8.4f</td>', $rd['z']);
            } elseif (abs((double)$arr[6]-$rd['z'])<0.001) {
                printf ('<td><font color="green">%8.4f</font></br>%8.4f</td>', $rd['z'],$arr[6]);
            } else {
                printf ('<td><font color="red">%8.4f</font></br>%8.4f</td>', $rd['z'],$arr[6]);
            }
            echo '</tr>';
        }
        $i++;
    }
    echo '</table>';
    echo $fout.' van '.$i.'<br>';
}
?>