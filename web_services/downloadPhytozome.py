#! /usr/bin/env python3
# AUTHOR: daniel.lang@helmholtz-muenchen.de

import xml.etree.ElementTree
import re
import sys
import requests
import getpass
import os.path

login	= sys.argv[1]
path	= sys.argv[2]
version	= sys.argv[3] if sys.argv[3] != None else 11

url_0	= 'https://signon.jgi.doe.gov/signon/create' 
url_1	= 'http://genome.jgi.doe.gov/ext-api/downloads/get-directory?organism=PhytozomeV%i'
url_2	= 'http://genome.jgi.doe.gov'
data 	= { 'login': login, 'password': getpass.getpass(prompt="your password:",stream=None) }

s = requests.session()
s.post(url_0, data)
r=s.get(url_1)
try: xml = xml.etree.ElementTree.fromstring(r.text)
except: error('Failed to parse returned XML')

for n in xml.findall("*//file[@filename]"):
	if n.attrib["filename"] == 'DataReleasePolicy.html': 
		continue
	if not os.path.exists("/".join([path,n.attrib["filename"]])):
		URL="/".join([url_2,n.attrib["url"]])
		print("Downloading %s" % URL)
		with open("/".join([path,n.attrib["filename"]]),"wb") as handle:
			req=s.get(URL,stream=True)
			for block in req.iter_content(1024):
				handle.write(block)

