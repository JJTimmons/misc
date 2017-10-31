proc loadAll {} {
	foreach n {0 1 2 3 4 5 6 7 8 9 10 11} {
		mol new ff${n}.pdb
	}
}

proc psfAll {} {
	# top is now the newly loaded mol
	mol new "protein.pdb"

	package require psfgen
	package require topotools
	resetpsf

	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_prot.rtf
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_lipid.rtf
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_na.rtf
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_carb.rtf
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/top_all36_cgenff.rtf
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/toppar_all36_carb_glycopeptide.str
	topology /usr/local/lib/vmd/plugins/noarch/tcl/readcharmmtop1.2/toppar_water_ions_namd.str
	topology /home/josh/Desktop/9_3_tubBendPMF/equils/0/parameters/toppar_all36_na_nad_ppi_gdp_gtp.str

	proc chainGen {chain sel} {
		# segment by selection into a new chain file
		set chainPDB ${chain}.pdb
		[atomselect top $sel] writepdb $chainPDB

		# segment on the newly created PDB
		segment $chain {pdb $chainPDB}
		coordpdb $chainPDB $chain
	}

	chainGen A "chain A"
	chainGen B "chain B"
	chainGen X "chain X"
	chainGen Y "chain Y"
	chainGen Z "chain Z"

	guesscoord
	writepdb "protein.pdb"
	writepsf "protein.psf"
}

# move a to b and add GTP and GDP and MG from a to b
# this should be called from the directory above the
# directories that hold a and b
proc nucMov {a b} {

	set chainAa [atomselect $a "chain A"]
	set chainAb [atomselect $b "chain A"]
	[atomselect $a "chain X or chain Z"] move [measure fit $chainAa $chainAb]

	set chainBa [atomselect $a "chain B"]
	set chainBb [atomselect $b "chain B"]
	[atomselect $a "chain Y"] move [measure fit $chainBa $chainBb]

	set sellist {}
	lappend sellist [atomselect $a "chain X or chain Y or chain Z"]
	lappend sellist [atomselect $b "protein"]

	set mol [::TopoTools::selections2mol $sellist]
	animate write psf protein.psf $mol
	animate write pdb protein.pdb $mol
}

# prepare by solvating and adding ions
# hardcoded center based on tubulin post molmov
# this should be run after the tubulin gets its
# GTP, GDP and MG
proc prepare {} {
	package require psfgen 1.5
	package provide solvate 1.7
	package require autoionize
	resetpsf

	solvate protein.psf protein.pdb -o protein -minmax {{-55 -77 -70} {95 73 80}}
	autoionize -psf protein.psf -pdb protein.pdb -neutralize -o protein

	set all [atomselect top all]
	set pro [atomselect top "protein and backbone"]
	$all set beta 0
	$pro set beta 1
	$all writepdb "protein_restraints.pdb"
}
