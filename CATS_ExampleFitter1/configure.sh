#!/bin/bash

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

if ! cmake .; then
	echo "Configuration of CATS ${red}failed${reset}"
	return 3
fi

echo "The configuration of CATS was ${green}successful${reset}"
#echo "  Configured using "$1
echo "  To proceed type: make"


return 0
