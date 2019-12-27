# Check if there is ring piercing in current top molecule in VMD
# by measuring all bond lengths longer than 2.0 Ang.
########################################
# Set name of PSF file.
set PSF "step5_charmm2namd.psf"

########################################
# Execute only if PSF exists
if {[file exists $PSF]} {
    set ID [molinfo top]
    set listpiercing {}
    set INP [open $PSF r]
    
    set line [gets $INP]
    while { ![string match -nocase "*NBOND: bonds*" $line] } {
        set line [gets $INP]
    }
    while { ![string match -nocase "" $line] } {
        set line [gets $INP]
        for {set i 1} {$i<=[expr [llength $line]/2]} {incr i} {
            set b [measure bond [list [expr [lindex $line [expr $i*2-2]] - 1]  [expr [lindex $line [expr $i*2-1]] - 1] ]]
           	if {$b>2.0 && $b<15} {
    	    lappend listpiercing [list [lindex $line [expr $i*2-2]]  [lindex $line [expr $i*2-1]]]
                puts "   Bond between serial [lindex $line [expr $i*2-2]]  [lindex $line [expr $i*2-1]] is [format %.3f $b]"
    
    	}
        }
    }
    close $INP
    if {![llength $listpiercing]} {puts "* No ring piercing found *"; return}
    # Check if there is a disulfide bond, which might have a bond longer than 2.0 A. Remove from the list in that case.
    set listfinal {}
    for {set i 0} {$i<[llength $listpiercing]} {incr i} {
        if {[string trim [[atomselect $ID "serial [lindex $listpiercing $i 0]"] get type]] eq "SM" && [string trim [[atomselect $ID "serial [lindex $listpiercing $i 1]"] get type]] eq "SM"} {
            puts "      Serial [lindex $listpiercing $i 0]  [lindex $listpiercing $i 0] (atom type \"SM\") is NOT a ring piercing, but rather a disulfide bond."
        } else {
            lappend listfinal [lindex $listpiercing $i]
        }
    }

} else {
    puts "PSF file does not exist!"
    return
}
# Procedure to add molecular representation
proc AddRepPiercing { SERIAL MOLID } {
    mol color Name
    mol representation VDW 0.8 12.0
    mol selection "serial $SERIAL"
    mol addrep $MOLID
    mol color ResName
    mol representation Licorice 0.1 12.0 12.0
    mol selection "same residue as(within 3 of (serial $SERIAL))"
    mol addrep $MOLID
}
# Add representation for piercing atoms
AddRepPiercing [join $listfinal] $ID

