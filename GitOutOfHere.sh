#!/bin/bash


# set it up so I cna push wiht an ssh agent (no user name or password)
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_rsa_git

# a script for the end of the day that pushes all my changes to the repo
git commit -a
git push -u origin master



