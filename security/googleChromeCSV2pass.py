#!/usr/bin/env python3

# shamelessly adapted by Daniel Lang <dandaman@couchbanditen.de> from source code of  
# Copyright 2015 David Francoeur <dfrancoeur04@gmail.com>
# This file is licensed under the GPLv2+. Please see COPYING for more information.

# Usage: ./googleChromeCSV2pass.py googleExport.csv [folder_name|chromePW]

import csv
import itertools
import sys
import os
import time
from subprocess import Popen, PIPE

def pass_import_entry(path, data):
	""" Import new password entry to password-store using pass insert command """
	proc = Popen(['pass', 'insert', '--multiline', path], stdin=PIPE, stdout=PIPE)
	proc.communicate(data.encode('utf8'))
	proc.wait()

def readFile(filename):
	""" Read the file and proccess each entry """
	with open(filename, 'rU') as csvIN:
		next(csvIN)
		outCSV=(line for line in csv.reader(csvIN, dialect='excel'))
		#for row in itertools.islice(outCSV, 0, 1):
		for row in outCSV:
			#print("Length: ", len(row), row) 
			prepareForInsertion(row)


def prepareForInsertion(row):
	""" prepare each CSV entry into an insertable string """
	try: 
		sys.argv[2]
	except IndexError:	
		keyFolder = "chromePW"
	else:
		keyFolder = escape(sys.argv[2])

	keyName = escape(row[0])
	username = row[2] if row[2] else "NO_USER"
	password = row[3]
	url = row[1]
	
	path = keyFolder+"/"+keyName+"/"+username
	data = password + "\n" if password else "\n"
	data = "%s%s: %s\n" % (data, "url:", url+"\n")
	data = "%s%s: %s\n" % (data, "imported on:", time.strftime("%c") +"\n")
	pass_import_entry(path,data) 
	print(path+" imported!")

def escape(strToEscape):
	""" escape the list """
	return strToEscape.replace(" ", "-").replace("&","and").replace("[","").replace("]","")
	

def main(argv):
	inputFile = sys.argv[1]
	print("File to read: " + inputFile)
	readFile(inputFile)
	
	yes = set(['yes','y', 'ye', ''])
	choice=input("Do you want to delete " + inputFile + "? [y|N]: ").lower()
	if choice in yes:
		os.remove(inputFile)
		print(inputFile+" deleted!\n")


main(sys.argv)
