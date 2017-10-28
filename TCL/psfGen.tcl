proc loadAll {} {
	foreach n {0 1 2 3 4 5 6 7 8 9 10 11} {
		mol new ff${n}.pdb
		mol addfile ff${n}.psf
	}
}

proc psfAll {} {
	package require psfgen
	package require topotools

	foreach n [molinfo list] {
		psfcontext reset
		resetpsf
		topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf
		set chainA A_$n.pdb
		[atomselect $n "chain A"] writepdb $chainA
		segment "A" {pdb $chainA}
		coordpdb $chainA "A"
		guesscoord

		set name [molinfo $n get "name"]
		set autoNameA ${name}_Aauto
		writepdb $autoNameA.pdb
		writepsf $autoNameA.psf

		########
		resetpsf
		topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf
		set chainB B_$n.pdb
		[atomselect $n "chain B"] writepdb $chainB
		segment "B" {pdb $chainB}
		coordpdb $chainB "B"
		guesscoord

		set autoNameB ${name}_Bauto
		writepdb $autoNameB.pdb
		writepsf $autoNameB.psf

		########
		set midlist {}
		set mol [mol new $autoNameA.pdb waitfor all]
		mol addfile $autoNameA.psf
		lappend midlist $mol
		set mol [mol new $autoNameB.pdb waitfor all]
		mol addfile $autoNameB.psf
		lappend midlist $mol
		set mol [::TopoTools::mergemols $midlist]
		animate write psf ${name}_auto.psf $mol
		animate write pdb ${name}_auto.pdb $mol
		lappend midList $mol

		foreach m $midList {
			mol delete $m
		}

		########
		file delete $chainA
		file delete $chainB
		file delete $autoNameA.pdb
		file delete $autoNameA.psf
		file delete $autoNameB.pdb
		file delete $autoNameB.psf
	}
}

# move a to b and add GTP and GDP and MG from a to b
proc nucMov {a b} {

	set neighborsA [atomselect $a "name CA and within 7 of (resname GTP or resname GDP)"]

	set query ""
	foreach r [$neighborsA get resid] {
		append query " or resid ${r}"
	}
	set query [string range $query 4 [expr [string length $query] - 1]]

	

}


namdenergy -sel [atomselect top "resid gt 434 and protein"] -vdw -ofile "out.dat" -par /home/josh/Desktop/10_11/parameters/par_all36_na.prm -par /home/josh/Desktop/10_11/parameters/par_all36_cgenff.prm -par /home/josh/Desktop/10_11/parameters/par_all36_carb.prm -par /home/josh/Desktop/10_11/parameters/par_all36_na.prm -par /home/josh/Desktop/10_11/parameters/par_all36_prot.prm -par /home/josh/Desktop/10_11/parameters/toppar_all36_carb_glycopeptide.str -par /home/josh/Desktop/10_11/parameters/toppar_all36_na_nad_ppi_gdp_gtp.str -par /home/josh/Desktop/10_11/parameters/toppar_water_ions_namd.str -par /home/josh/Desktop/10_11/parameters/par_all36_lipid.prm