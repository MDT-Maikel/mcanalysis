#! /usr/bin/env python

# python modules
import sys
import numpy

# mcanalysis python module
import mca


# program arguments
if len(sys.argv) < 2:
	print "Usage: " + sys.argv[0] + "<event dir> <analysis program>"
	sys.exit(1)
event_dir = sys.argv[1]
analysis_exe = "/home/maikel/programs/mcanalysis/build/analyses/dijet_cms_exo-11-016" #sys.argv[2]

outputfile = open("test.txt", "w")

masses_m = numpy.arange(250.0, 1400.0, 100.0)
for mass_m in masses_m:
	event_file = event_dir + "/mass_" + str(mass_m) + ".lhco.gz"
	acceptance = mca.run_analysis(analysis_exe, event_file)
	xsec = mca.read_mg5_xsec(event_dir + "/mass_" + str(mass_m) + ".txt")
	outputstring = "%.3e %.3e %.3e\n" % (mass_m, xsec, float(acceptance))
	outputfile.write(outputstring)

