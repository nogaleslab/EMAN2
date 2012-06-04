#!/usr/bin/env python
#test added by mike
import os
from global_def import *
from optparse import OptionParser
import sys

from EMAN2  import *

def main():
	arglist = []
        for arg in sys.argv:
		arglist.append( arg )
	        progname = os.path.basename(arglist[0])
	        usage = progname + " raw_images apix <CTFparfile>"
	        parser = OptionParser(usage,version="1.0")
	        (options, args) = parser.parse_args(arglist[1:])
        if len(args) < 2:
                print "usage: " + usage
                print "Please run '" + progname + " -h' for detailed options"
	else:
		a = EMData()
		fname = args[0]
		apix = float(args[1])
		ctfpar = None
		if len(args) == 3:
			ctfpar = args[2]
			if not os.path.isfile(ctfpar):
				print "Error: ctf file '%s' is not found"%ctfpar
				sys.exit()
			# get ctf values from parfile
			from utilities import generate_ctf
			pfile = open(ctfpar).readlines()
			parts=[0]*len(pfile)
			for p in pfile:
				nfo = p.strip().split()
				pnum=int(float(nfo[0]))-1
				df1=float(nfo[1])
				df2=float(nfo[2])
				astig = float(nfo[3])
				kv = int(float(nfo[4]))
				cs = float(nfo[5])
				ampc = float(nfo[6])*100
				
				if len(nfo) == 9:
					hnum = int(float(nfo[7]))
					angle = int(float(nfo[8]))

				"""
				# For FREALIGN
				if p[0]=="C":
					continue
				(pnum,psi,theta,phi,shx,shy,mag,hnum,ctf1,ctf2,angle,pres,delta)=p.strip().split()
				pnum=int(pnum)-1
				kv=300
				cs=2.7
				ampc=5
				"""
				df=(float(df1)+float(df2))/2
				parts[pnum]=generate_ctf([df,cs,kv,apix,0.0,ampc])

		imn = EMUtil.get_image_count(fname) 
		print "Generating 'start.hdf' with %i particles"%imn
		for i in xrange(imn):
			a.read_image(fname,i)
			a.set_attr_dict({'active':1})
			if ctfpar is not None:
				a.set_attr("ctf",parts[i])
			t2 = Transform({"type":"spider","phi":0,"theta":0,"psi":0})
		        a.set_attr("xform.projection", t2)
			a.set_attr("apix_x",apix )
			a.write_image("start.hdf",i)
			print "%3i%% complete\t\r"%(int(float(i)/imn*100)),
		print "100% complete\t"

if __name__ == "__main__":
        main()

