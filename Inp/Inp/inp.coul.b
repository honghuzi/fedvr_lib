$route

$end
          
$nonstd
11//00;
20//01;
$end

$title
DVR Code
$end

$dvr_input
 number-of-space-variables=1 space-variable-1=x
 kinetic-energy-type=packed angular-momentum=0
$end

$h0(x)
 automate number-of-major-blocks=3 
 left-fixed-point right-fixed-point reuse-space-data
 drop-left-function drop-right-function print=all
$end

 $block-1
  default-order=(50) number-of-subregions=1
  left-boundary=1.d-10 right-boundary=20.d0
 $end

 $block-2
  default-order=(15) number-of-subregions=1
  left-boundary=20.d0 right-boundary=100.d0
 $end

 $block-3
  default-order=(10) number-of-subregions=1
  left-boundary=100.d0 right-boundary=200.d0
 $end

$v_reg_1(x)
 use-atomic-units potential=coulomb 
$end
