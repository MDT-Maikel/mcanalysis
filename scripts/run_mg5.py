#! /usr/bin/env python

# python modules
import sys
import numpy

# mcanalysis python module
import mca

# program arguments
if len(sys.argv) < 4:
	print "Usage: " + sys.argv[0] + "<process dir> <event dir> <output file>"
	sys.exit(1)
proc_dir = sys.argv[1]
event_dir = sys.argv[2]
output_file = sys.argv[3]
param_card_file = proc_dir + "/Cards/param_card.dat"

outputfile = open(output_file, "w")

# Grid loops
masses = numpy.arange(200.0, 4000.0, 200.0)
for mass in masses:
	mca.change_mg5_par(param_card_file, "MM", mass)
	mca.update_width_mg5(proc_dir)
	mca.run_mg5(proc_dir, event_dir, "mass_" + str(mass), True)
	xsec = mca.read_mg5_xsec(event_dir + "/mass_" + str(mass) + ".txt")
	outputstring = "%.3e %.3e\n" % (mass, xsec)
	outputfile.write(outputstring)


#print read_mg5_xsec("/home/maikel/programs/mcanalysis/scripts/test/test.txt")

#print read_mg5_par(proc_dir + "/Cards/param_card.dat", "yLu3x3")

#change_mg5_par(proc_dir + "/Cards/param_card.dat", "yLu3x3", 0.123454)


	
