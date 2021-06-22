#!/bin/bash

find . \( -type d -and \( -name "log" -or -name "logs" \) \) -prune > log_dir_locations.txt





find . \( -type d -and \( -name "unitigging" -or -name "canu-logs" -or -name "unitigging.html.files" -or -name "trimming.html.files" -or -name "correction.html.files" -or -name "canu-scripts"  \) \) -prune -print

						# too generic to delete 
						#-name "correction" -or \
						#-name "trimming" -or \
