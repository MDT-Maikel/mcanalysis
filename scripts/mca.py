#! /usr/bin/env python
##########################################################
## This module contains useful functions for generating ##
## Monte Carlo data using MadGraph and analyzing it     ##
## using the functionality of this analysis suite       ##
##########################################################

# python modules
import sys
import os
import shutil
import random
import subprocess
import re


########################################################
## runs an instance of madgraph and stores the events ##
########################################################
def run_mg5(proc_dir, event_dir, run_name, delphes = False):
	# determine a unique process name
	proc_name = "run_" + str(random.randrange(0, 2**28))

	# check if madevent exists in the process dir
	if not os.path.isfile(proc_dir + "/bin/madevent"):
		print "Error: madevent is not available in process dir (" + proc_dir + ")."
		sys.exit(1)

	# write run details into a file
	mg5_run_details = "generate_events " + proc_name + "\npythia=OFF\ndelphes=OFF\n0"
	if delphes:
		mg5_run_details = "generate_events " + proc_name + "\npythia=ON\ndelphes=ON\n0"
	with open(proc_dir + "/mg5_script", 'w') as mg5_run_file:
		mg5_run_file.write(mg5_run_details)

	# run madgraph with the run details
	os.system(proc_dir + "/bin/madevent < " + proc_dir + "/mg5_script")
	os.remove(proc_dir + "/mg5_script")

	# create the event dir
	if not os.path.isdir(event_dir):
		os.makedirs(event_dir)

	# store the results in the event dir
	if os.path.isfile(proc_dir + "/Events/" + proc_name + "/" + proc_name + "_tag_1_banner.txt"):
		shutil.move(proc_dir + "/Events/" + proc_name + "/" + proc_name + "_tag_1_banner.txt", event_dir + "/" + run_name + ".txt")
	if os.path.isfile(proc_dir + "/Events/" + proc_name + "/unweighted_events.lhe.gz"):
		shutil.move(proc_dir + "/Events/" + proc_name + "/unweighted_events.lhe.gz", event_dir + "/" + run_name + ".lhe.gz")
	if os.path.isfile(proc_dir + "/Events/" + proc_name + "/tag_1_delphes_events.lhco.gz"):
		shutil.move(proc_dir + "/Events/" + proc_name + "/tag_1_delphes_events.lhco.gz", event_dir + "/" + run_name + ".lhco.gz")
	if os.path.isfile(proc_dir + "/HTML/" + proc_name + "/results.html"):
		shutil.move(proc_dir + "/HTML/" + proc_name + "/results.html", event_dir + "/" + run_name + ".html")
	
	# return nothing
	return


#################################################
## updates the width in a madgraph process dir ##
#################################################
def update_width_mg5(proc_dir):
	# check if madevent exists in the process dir
	if not os.path.isfile(proc_dir + "/bin/madevent"):
		print "Error: madevent is not available in process dir (" + proc_dir + ")."
		sys.exit(1)

	# write run details into a file
	mg5_run_width_details = "compute_widths all --body_decay=2"
	with open(proc_dir + "/mg5_script_width", 'w') as mg5_run_width_file:
		mg5_run_width_file.write(mg5_run_width_details)

	# run madgraph with the run details
	os.system(proc_dir + "/bin/madevent < " + proc_dir + "/mg5_script_width")
	os.remove(proc_dir + "/mg5_script_width")

	# return nothing
	return


###########################################
## reads a madgraph cross section result ##
###########################################
def read_mg5_xsec(banner_file):
	# check if file exists
	if not os.path.isfile(banner_file):
		print "Error: banner file for reading cross section not found (" + banner_file + ")."
		return 0
	# read the banner file as a string	
	with open(banner_file, "r") as read_banner:
		banner = read_banner.read()
	# find the cross section
	xsec_pos_begin = banner.find("<MGGenerationInfo>") + len("<MGGenerationInfo>") + 1
	xsec_pos_end = banner.find("</MGGenerationInfo>") - 1
	banner = banner[xsec_pos_begin:xsec_pos_end]
	banner = banner.split()	
	banner = banner.pop()
	# convert and return the cross section
	return float(banner)


###################################################
## change a madgraph parameter in the param_card ##
###################################################
def change_mg5_par(param_card_file, param_id, value):
	# check if file exists
	if not os.path.isfile(param_card_file):
		print "Error: param card for writing parameter not found (" + param_card_file + ")."
		return
	# read the param card as a string	
	with open(param_card_file, "r") as read_card:
		param_card = read_card.read()
	# find the line for this parameter
	param_identifier = "# " + param_id + " "
	param_line = [line for line in param_card.split("\n") if param_identifier in line].pop()
	param_line_pos = param_line.find(param_identifier)
	param_line_val = param_line[0:param_line_pos]
	param_line_val = param_line_val.split().pop()
	param_line_new = param_line.replace(param_line_val, str(value))
	new_param_card = param_card.replace(param_line, param_line_new)
	# write the new param card
	with open(param_card_file, "w") as write_card:
		write_card.write(new_param_card)
		write_card.truncate()
	# return nothing
	return


###################################################
## read a madgraph parameter from the param_card ##
###################################################
def read_mg5_par(param_card_file, param_id):
	# check if file exists
	if not os.path.isfile(param_card_file):
		print "Error: param card for reading parameter not found (" + param_card_file + ")."
		return 0
	# read the param card as a string	
	with open(param_card_file, "r") as read_card:
		param_card = read_card.read()
	# find the parameter
	param_pos = param_card.find("# " + param_id)
	param_card = param_card[max(0, param_pos - 40):param_pos]
	param_card = param_card.split().pop()
	# return the parameter value
	return float(param_card)


#################################################
## change a madgraph parameter in the run_card ##
#################################################
def change_mg5_run_par(run_card_file, param_id, value):
	value = str(value)
	# check if file exists
	if not os.path.isfile(run_card_file):
		print "Error: run card for writing parameter not found (" + run_card_file + ")."
		return
	# read the run card as a string	
	with open(run_card_file, "r") as read_card:
		run_card = read_card.read()
	param_identifier = "[ \t]*([a-zA-Z0-9.]+)[ \t]+=[ \t]+" + param_id + "[ \t]+"
	# check if identifier could be found
	if not re.search(param_identifier, run_card):
		print "Error: param identifier could not be found in run card (" + run_card_file + ")."
		return
	# substitute the new value
	new_run_card = re.sub(param_identifier, lambda m: re.sub(m.group(1), value, m.group(0)), run_card)
	# write the new param card
	with open(run_card_file, "w") as write_card:
		write_card.write(new_run_card)
		write_card.truncate()
	# return nothing
	return


############################
## run an analysis script ##
############################
def run_analysis(analysis_program, event_file):
	return_value = subprocess.check_output([analysis_program + " " + event_file], shell = True)
	acceptances = [float(x) for x in return_value.split()]
	return acceptances



