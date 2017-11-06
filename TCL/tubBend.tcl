
# {param} m1 - Number [id of the reference molecule of 1JFF, must have X and Y chains]
# {param} m2 - Number [id of the molecule it's being compared against]
proc tubBend {m1 m2} {

	# set 1JFF (ref), chain Z, chain Y
	set m1a A
	set m1b B

	# set test mol, chain A, chain B
	set m2a A
	set m2b B
	# set m2a C
	# set m2b D
	# set m2a X
	# set m2b Y


	proc tubBendCalc {m1 m1a m1b m2 m2a m2b} {
		#set all H7 selection keywords
		set a1 "chain ${m1a} and resid 222 to 244 and backbone"
		set b1 "chain ${m1b} and resid 222 to 244 and backbone"
		set a2 "chain ${m2a} and resid 222 to 244 and backbone"
		set b2 "chain ${m2b} and resid 222 to 244 and backbone"

		#set the selections
		set 1H7a [atomselect $m1 $a1]
		set 1H7b [atomselect $m1 $b1]
		set 2H7a [atomselect $m2 $a2]
		set 2H7b [atomselect $m2 $b2]

		#find fit between the two alpha H7 helixes
		set 1fit [measure fit $1H7a $2H7a]

		#realign full molecules using fit
		[atomselect $m1 all] move $1fit

		#find helix centers
		set A1c [measure center $1H7a weight mass]
		set B1c [measure center $1H7b weight mass]
		set A2c [measure center $2H7a weight mass]
		set B2c [measure center $2H7b weight mass]

		#create vectors between each mol's H7's
		set v1 [vecsub $B1c $A1c]
		set v2 [vecsub $B2c $A2c]

		#divide dot product by multiplied magnitudes
		set dm [expr [vecdot $v1 $v2] / [expr [veclength $v1] * [veclength $v2]]]

		#arccos
		return [expr 57.2957795 * [expr acos($dm)]]
	}

	# open file
	set nm [molinfo $m2 get name]
	set output [open "${nm}_bend.dat" w]

	# bend calculation loop
	set n [molinfo $m2 get numframes]
	for {set i 0} {$i < $n} {incr i} {
		molinfo $m2 set frame $i
		set bend [tubBendCalc $m1 $m1a $m1b $m2 $m2a $m2b]
		puts $bend
		# puts $output "$bend"
	}
	# puts "output file: ${nm}_bend.dat"
	close $output
}
