#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

if ! cmake .; then
	echo "Configuration ${red}failed${reset}"
	return 3
fi
echo "The configuration was ${green}successful${reset}"
echo "  To proceed type: make"

return 0
