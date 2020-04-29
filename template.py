#!/usr/bin/env python
import argparse
import os 
import sys

# global var for inputs
args=None 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("-s", "--string", help="string option")
	parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()

