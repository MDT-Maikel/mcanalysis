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
analysis_exe = sys.argv[2]

outputfile = open("test.txt", "w")

masses_m = numpy.arange(200.0, 1400.0, 100.0)
for mass_m in masses_m:
	event_file = event_dir + "/mass_" + str(mass_m) + ".lhco.gz"
	acceptance = mca.run_analysis(analysis_exe, event_file)
	xsec = mca.read_mg5_xsec(event_dir + "/mass_" + str(mass_m) + ".txt")
	result_list = [mass_m, xsec] + acceptance
	outputstring = "";
	for i in result_list:
		outputstring = outputstring + ("%.3e " % i)
	outputstring = outputstring + "\n";
	outputfile.write(outputstring)

