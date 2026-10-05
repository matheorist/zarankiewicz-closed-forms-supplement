# =====================================================================
#  Supplementary numerical verification for
#  "Closed-form expressions for a Zarankiewicz linear program"
#  (Y. Zhang, G. Chen)
#
#  Float-based checks of the paper's quantitative claims.  (The exact
#  rational simplex / CAS identities are a separate engine; this file
#  reproduces the numbers a reader can spot-check.)  Run:  pwsh verification.ps1
# =====================================================================

$cbrt2 = [math]::Pow(2, 1.0/3)   # 2^(1/3) = 1.259921...

# ---------- closed forms ----------
function Fq($q,$a,$b){
  $U=4.0*$a*$a*$a-4*$a*$a*$q-6*$a*$q*$q*$q+9*$a*$q*$q-6*$q*$q*$q
  $V=2.0*$a*$a*$a+5*$a*$a*$q-4*$a*$a-12*$a*$q*$q*$q+12*$a*$q*$q+26*$a*$q-30*$q*$q*$q-6*$q*$q
  $W=4.0*$a*$a*$a*$a-12*$a*$a*$q*$q*$q+17*$a*$a*$q*$q-2*$a*$q*$q*$q*$q-15*$a*$q*$q*$q+4*$q*$q*$q*$q
  (2*(2*$a-$q)*($a+$q)*$U*$b+1.0*$a*$q*$q*$q*($a-1)*$V)/(3*$q*$W)
}
function Bminus($q,$a){ [double]$q*$q*$q*($a-1)*($a+4)/(2*($a+$q)*(2*$a+$q)) }
function Bplus($q,$a){ [double]$a*$q*$q*($a-1)*(2*$a*$q+6*$a-9*$q*$q-$q)/(2*(2*$a-$q)*(2*$a-3*$q)*($a-$q)) }
function C2([double]$x){ $x*($x-1)/2.0 }
function Fqt([int]$q,[int]$a,[double]$b){
  $c=2*$a/$q
  if($c -ne [math]::Floor($c)){ throw "Fqt requires q | 2a" }
  $c=[int]$c
  $t=[math]::Floor(($c-$q+2)/($q+1))
  if($t -lt 2){ throw "point is outside the t>=2 tiles" }
  $u=$c-$t
  $r=$c-(($q+1)*$t+$q-3)
  if(($r -lt 1)-or($r -gt $q+1)){ throw "point is outside the migrated tile" }
  $C=C2 $a
  $ncm1=4*$C/(C2 ($c+1))
  $ncm2=($q-4)*$C/(C2 ($c+2))
  $nu=$r*$C/(C2 $u)
  $nup1=($q+1-$r)*$C/(C2 ($u+1))
  $bm=$ncm1+$ncm2
  $bp=$nu+$nup1
  $Zm=($c+1)*$ncm1+($c+2)*$ncm2
  $Zp=$u*$nu+($u+1)*$nup1
  $lambda=($b-$bm)/($bp-$bm)
  (1-$lambda)*$Zm+$lambda*$Zp
}

# ---------- Roman's bound R(V,E): 2-constraint LP, optimum on <=2 sizes ----------
function Comb3([int]$i){ if($i -lt 3){0.0} else {[double]$i*($i-1)*($i-2)/6.0} }
function Roman([int]$V,[int]$E){
  $budget=2.0*(Comb3 $V); $best=0.0
  for($i=0;$i -le $V;$i++){
    $ci=Comb3 $i
    if(($ci*$E) -le $budget){ $val=[double]$i*$E; if($val -gt $best){$best=$val} }
    for($j=$i+1;$j -le $V;$j++){
      $cj=Comb3 $j
      if($cj -ne $ci){ $nj=($budget-$ci*$E)/($cj-$ci); $ni=[double]$E-$nj
        if(($ni -ge 0.0)-and($nj -ge 0.0)){ $val=$i*$ni+$j*$nj; if($val -gt $best){$best=$val} } }
    }
  }
  ,$best
}

Write-Output "=== (1) Paper comparison table: F_q, floor F_q, floor Roman, integer gain ==="
foreach($r in @(5,40,40),@(5,45,40),@(6,60,60),@(6,66,66),@(7,70,84),@(7,70,86)){
  $q=$r[0];$a=$r[1];$b=$r[2]; $F=Fq $q $a $b; $R=Roman $a $b
  "  q={0} a={1} b={2}:  F_q={3:N2}  floorF={4}  floorR={5}  gain={6}" -f `
     $q,$a,$b,$F,[math]::Floor($F),[math]::Floor($R),([math]::Floor($R)-[math]::Floor($F))
}

Write-Output "`n=== (2) Dual-certificate identity:  dual objective == F_q ==="
function DualObj($q,$a,$b){
  $W=4.0*$a*$a*$a*$a-12*$a*$a*$q*$q*$q+17*$a*$a*$q*$q-2*$a*$q*$q*$q*$q-15*$a*$q*$q*$q+4*$q*$q*$q*$q
  $U=4.0*$a*$a*$a-4*$a*$a*$q-6*$a*$q*$q*$q+9*$a*$q*$q-6*$q*$q*$q
  $Nq=4.0*$a*$a*$a-24*$a*$a*$q+12*$a*$q*$q*$q-7*$a*$q*$q+18*$q*$q*$q*$q+3*$q*$q*$q
  $y0=2*(2*$a-$q)*($a+$q)*$U/(3*$q*$W)
  $ys=$q*$q*$Nq/((2*$a+$q)*$W)
  $yc1=2*$q*$q*$q*(2*$a-3*$q*$q+$q)*(10.0*$a*$a-5*$a*$q-6*$q*$q*$q+7*$q*$q)/(($a-$q)*(2*$a+$q)*$W)
  $yc=2*$q*$q*($a-$q*$q+$q)*(12.0*$a*$a*$a-32*$a*$a*$q+11*$a*$q*$q+18*$q*$q*$q*$q-9*$q*$q*$q)/(($a-$q)*(2*$a+$q)*$W)
  function C($n,$k){$r=1.0;for($t=0;$t -lt $k;$t++){$r*=($n-$t)};for($t=1;$t -le $k;$t++){$r/=$t};$r}
  $b*$y0+2*(C $a 3)*$ys+$q*(C $a 2)*($yc1+$yc)
}
foreach($r in @(5,40,40),@(6,60,60),@(7,70,84)){
  $q=$r[0];$a=$r[1];$b=$r[2]; $d=DualObj $q $a $b; $f=Fq $q $a $b
  "  q={0} a={1} b={2}:  dual={3:N5}  F_q={4:N5}  equal={5}" -f $q,$a,$b,$d,$f,([math]::Abs($d-$f) -lt 1e-6)
}

Write-Output "`n=== (3) Diagonal barrier:  F_(q,t)(n,n)/n^(5/3) -> 2^(1/3),  F_(q,t)/Roman -> 1 ==="
foreach($r in @(5,40),@(7,84),@(8,120),@(10,240),@(12,420),@(15,825)){
  $q=$r[0];$n=$r[1]; $F=Fqt $q $n $n; $R=Roman $n $n
  "  q={0} n={1}:  F_(q,t)/n^(5/3)={2:N5}  F_(q,t)/R={3:N5}" -f $q,$n,($F/[math]::Pow($n,5.0/3)),($F/$R)
}
"  (target: F_(q,t)/n^(5/3) -> 2^(1/3) = {0:N5};  F_(q,t)/R -> 1)" -f $cbrt2

Write-Output "`n=== (4) Brown elliptic-quadric construction (K33-free lower bound) ==="
function BlockSize($p,$d){ $c=0
  for($x=0;$x -lt $p;$x++){for($y=0;$y -lt $p;$y++){for($z=0;$z -lt $p;$z++){
    if((($x*$x+$y*$y+$z*$z)%$p) -eq $d){$c++} }}}; $c }
function IsRes($d,$p){ for($x=1;$x -lt $p;$x++){ if(($x*$x)%$p -eq $d){return $true} }; $false }
foreach($p in 5,7,11,13,17,19,23,29){
  # elliptic class: pick d giving p^2-p points
  $d=1; while((BlockSize $p $d) -ne ($p*$p-$p)){ $d++ ; if($d -ge $p){break} }
  $bs=BlockSize $p $d; $c0=$bs/($p*$p*1.0)
  "  p={0,3}  elliptic d={1,3}  |block|={2,5} (=p^2-p)  c0=z/N^(5/3)={3:N4} (=1-1/p)" -f $p,$d,$bs,$c0
}
Write-Output "  => z(N,N;3,3) >= (1-1/p) N^(5/3) -> N^(5/3);  upper bound floor F_(q,t) ~ 2^(1/3) N^(5/3)."
Write-Output "  => diagonal gap = 2^(1/3) ~ 1.260."

Write-Output "`n=== (5) Gap = block-size ratio:  LP block c=2a/q ~ 2^(1/3) a^(2/3) vs Brown ~ a^(2/3) ==="
foreach($r in @(84,7),@(240,10),@(825,15),@(5460,28)){
  $a=$r[0];$q=$r[1]; $c=2.0*$a/$q
  "  a={0,5} q={1,3}:  c/a^(2/3)={2:N4}  (-> 2^(1/3)={3:N4})" -f $a,$q,($c/[math]::Pow($a,2.0/3)),$cbrt2
}

Write-Output "`n=== (6) Exact integer feasibility checks for the two-edge examples ==="
python scripts/near_endpoint_gap_verify.py
if($LASTEXITCODE -ne 0){ exit $LASTEXITCODE }
