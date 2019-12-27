### SCRIPT ###
# Calculation of static dielectric constant for a solvent.
# The system should contain only the solvent molecules.
########     USER INPUT      #################
set stride 10  ;# skip frames in trajectory
set temp  298  ;# set the temperature
set start 200  ;# skip initial frames

##############################################
# Set initial variables:
set M []
set listMX []
set listMY []
set listMZ []
set listMX2 []
set listMY2 []
set listMZ2 []
set V []
set conv2D  0.20819434   ;# 1D = 0.20819434 e*Ang
set conv2CM 3.33564      ;# modified for efficiency, 1D = 3.33564x10^-30 C·m 
set kb      1.38064852   ;# modified for efficiency, 1.38064852x10^-23
set eps0    8.8541878176 ;# modified for efficiency, 8.8541878176x10-12
set frames []
# Loop over trajectory:
set numfr [molinfo top get numframes]
for {set f $start} {$f<$numfr} {set f [expr $f + $stride ]} {
    puts "FRAME $f"
    lappend frames $f
    animate goto $f
    set dipvec []
    set dipnorm []
    #set dipvecCOM []
    set sum {0. 0. 0.}
    # Loop over all solvent molecules:
    # Select all molecules in the frame
    set allsolv [atomselect top "all" frame $f]
    # Get the residue ID for each molecule
    set restmp [$allsolv get resid]
    set res  [lsort -unique -integer $restmp]
    # Get the number of residues
    set numres [llength $res]
    # Enter the loop:
    foreach n $res {
        set solv [atomselect top "resid ${n}" frame $f]
        set co [$solv get {x y z}]
        set ch [$solv get charge]
        set COM [measure center $solv weight mass]
        set NA [$solv num]
        set vectmp {0. 0. 0.}
        # Multiply each coordinate of each atom for the partial charge.
        set lis []
        for {set i 0} {$i<$NA} {incr i} {
            lappend lis [vecscale [lindex $ch $i] [lindex $co $i]]
            set vectmp [vecadd $vectmp [lindex $lis $i]]
        }  
        # Calculate dipole moment for the current solvent molecule
        lappend dipnorm [expr [veclength $vectmp]/$conv2D]
        #lappend dipvecCOM [vecadd [lindex $dipvec end] $COM];# Vector with respect to COM, for drawing
        # Calculate the vector of the total dipole moment
        set sum [vecadd $sum $vectmp]
    }
    lappend listMX  [lindex $sum 0]
    lappend listMY  [lindex $sum 1]
    lappend listMZ  [lindex $sum 2]
    # Calculate the magnitude of the total dipole moment
    lappend M [expr [veclength $sum]/$conv2D]
    puts "Average dipole moment = [format %.3f [expr [vecsum $dipnorm]/$numres]] D"
    puts "Total dipole moment   = [format %.3f [lindex $M end]] D"
    # Find volume of current frame in Ang^3
    set tmp [pbc get]
    lappend V [expr [lindex $tmp 0 0]*[lindex $tmp 0 1]*[lindex $tmp 0 2]]
}
# Number of frames used according to the stride
set NP [llength $M]
# Average volume.
set avV [expr ([vecsum $V])/$NP]
# Calculate fluctuations along the trajectory:
set MX [expr ([vecsum $listMX]/$NP)/$conv2D]
set MY [expr ([vecsum $listMY]/$NP)/$conv2D]
set MZ [expr ([vecsum $listMZ]/$NP)/$conv2D]
set avM1 [expr ($MX**2 +  $MY**2 +  $MZ**2)] ;#  <M>^2

set MX2 0
set MY2 0
set MZ2 0
for {set t 0} {$t < $NP} {incr t} {
   set MX2 [expr $MX2 + ([lindex $listMX $t]/$conv2D)**2 ]
   set MY2 [expr $MY2 + ([lindex $listMY $t]/$conv2D)**2 ]
   set MZ2 [expr $MZ2 + ([lindex $listMZ $t]/$conv2D)**2 ]
}
set MX2 [expr $MX2/$NP]
set MY2 [expr $MY2/$NP]
set MZ2 [expr $MZ2/$NP]
set avM2 [expr $MX2 +  $MY2 +  $MZ2] ;#  <M^2>

puts " "
puts "*******************************************"
puts "Using $NP frames for averages..."
puts "Squared average of M (<M>^2) = [format %.2f $avM1 ] D^2"
puts "Average of squared M (<M^2>) = [format %.2f $avM2 ] D^2"
puts "  <M^2>  -  <M>^2     = [format %.2f [expr $avM2 - $avM1]] D^2"
puts " "
# Static dielectric constant:
set eps [expr 1 + (10**5)*((($avM2 - $avM1)*($conv2CM**2))/(3*$eps0*$avV*$kb*$temp))]
puts "Static dielectric constant = [format %.3f $eps]"
puts " "
puts "*******************************************"

return
