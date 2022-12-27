#!/bin/bash

#downloads java packages and sets them up.

download_dir='/home/kenmn/bioinfo_tools/java_tools/'

package_links='java_pack.txt'

while read p; do

	wget --directory-prefix=$download_dir $p 

done < $package_links

unzip '$download_dir'*
