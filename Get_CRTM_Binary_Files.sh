#https://bin.ssec.wisc.edu/pub/s4/CRTM/fix_REL-3.0.0_20230802.tgz
foldername="fix_REL-3.0.0_20230802"

filename="${foldername}.tgz"
if test -f "$filename"; then
    if [ -d "fix/" ]; then #fix directory exists
        echo "fix/ already exists, doing nothing."
    else
        #untar the file and move directory to fix
				tar -zxvf $filename
				cd $foldername/
				mv fix ..
				cd ..
				rmdir $foldername
				echo "fix/ directory created from existing $filename file."
    fi 
else
    #download, untar, move
		echo "Downloading $filename, please wait about 5 minutes (4 GB tar file)"
	  wget  https://bin.ssec.wisc.edu/pub/s4/CRTM/$filename # CRTM binary files, add "-q" to suppress output. 
		
	  tar -zxvf $filename
		cd $foldername/
		mv fix ..
		cd ..
		rmdir $foldername
		echo "fix/ directory created from downloaded $filename."
fi
echo "Completed."
