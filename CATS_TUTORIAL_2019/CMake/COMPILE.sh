#!/bin/bash
pathpat="(/[^/]*)+:[0-9]+"
ccred=$(echo -e "\033[0;41m")
ccyellow=$(echo -e "\033[0;43m")
ccend=$(echo -e "\033[0m")
make -j4 "$@" 2>&1 | sed -E -e "/[Ee]rror[: ]/ s%$pathpat%$ccred&$ccend%g" -e "/[Ww]arning[: ]/ s%$pathpat%$ccyellow&$ccend%g"


