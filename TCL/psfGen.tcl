proc loadAll {} {
	foreach n {0 1 2 3 4 5 6 7 8 9 10 11} {
		mol new ff${n}.pdb
	}
}

proc psfAll {n} {
	package require psfgen
	package require topotools
	resetpsf

	# top is now the newly loaded mol
	mol new $n

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
	writepdb protein_3.pdb
	writepsf protein_3.psf
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
	animate write psf protein_merge.psf $mol
	animate write pdb protein_merge.pdb $mol
}

# prepare by solvating and adding ions
# hardcoded center based on tubulin post molmov
# this should be run after the tubulin gets its
# GTP, GDP and MG
proc prepare {} {
	package require psfgen 1.5
	package provide solvate 1.7
	package require autoionize

	solvate protein.psf protein.pdb -o protein -minmax {{-55 -77 -70} {95 73 80}}
	autoionize -psf protein.psf -pdb protein.pdb -neutralize -o protein

	set all [atomselect top all]
	set pro [atomselect top "protein and backbone"]
	$all set beta 0
	$pro set beta 1
	$all writepdb "protein_restraints.pdb"
}

# for looping thru and building each
foreach n {4 5 6 7} {
	cd $n
	if {$n != 10} {
		# load the last images
		set m [mol new "../[expr {$n - 1}]/protein.psf"]
		mol addfile "../[expr {$n - 1}]/npt.restart.coor"
		set o [mol new protein.pdb]

		# move the nucleotides and MG
		nucMov $m $o

		# PSF generate
		psfAll

		#prepare for MD
		prepare
	}

	foreach n [molinfo list] {
		mol delete $n
	}

	namd2 $sets $min
	namd2 $sets $nvt
	namd2 $sets $npt
	cd ".."
}

# look thru each to check for current bend
foreach n [molinfo list] {
	mol delete $n
}
set m [mol new 0/protein.psf]
mol addfile 0/npt.restart.coor
foreach n {1 2 3 4 5 6 7 8 9 10} {
	set o [mol new ${n}/protein.psf]
	mol addfile ${n}/npt.restart.coor
	
	set a [atomselect $m "name CA and (resid 206 to 215 or resid 224 to 242 or resid 252 to 259 or resid 290 to 300 or resid 325 to 335 or resid 269 to 272 or resid 312 to 320 or resid 352 to 356 or resid 377 to 381)"]
	set b [atomselect $o "name CA and (resid 206 to 215 or resid 224 to 242 or resid 252 to 259 or resid 290 to 300 or resid 325 to 335 or resid 269 to 272 or resid 312 to 320 or resid 352 to 356 or resid 377 to 381)"]

	puts [measure rmsd $a $b]
}