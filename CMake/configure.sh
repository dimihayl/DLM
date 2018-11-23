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
echo "                   make install"

CATSSYS=$(dirname $(pwd))
echo ""
echo "To have access to the CATS libraries you will need to run (and/or add in your .bashrc):"
echo '  export LD_LIBRARY_PATH=$(<PATH_TO_CATS>/bin/cats-config --libdir)${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}'

return 0
