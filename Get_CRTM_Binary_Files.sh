#https://bin.ssec.wisc.edu/pub/s4/CRTM/fix_REL-3.1.0.1_01252024.tgz  (use this for jedi and stand-alone, some files have changed).

# This script is used to manually download the tarball of binary and netcdf coefficient files.
# The same files also download automatically during the cmake step, so you don't have to actually run this manually. 

foldername="fix_REL-3.1.0.1_01252024"
filename="${foldername}.tgz"
echo "$filename"
break

if test -f "$filename"; then
    if [ -d "fix/" ]; then #fix directory exists
        echo "fix/ already exists, doing nothing."
    else
        #untar the file and move directory to fix
				tar -zxvf $filename
				mv $foldername/fix .
				rm -rf $foldername
				echo "fix/ directory created from existing $filename file."
    fi 
else
    #download, untar, move
		echo "Downloading $filename, please wait about 5 minutes (4 GB tar file)"
	  wget  https://bin.ssec.wisc.edu/pub/s4/CRTM/$filename # CRTM binary files, add "-q" to suppress output. 
				
    #untar the file and move directory to fix
    tar -zxvf $filename
    mkdir fix
    mv crtm/$foldername/* fix/.
    rm -rf crtm/$foldername
  	echo "fix/ directory created from downloaded $filename."
fi
echo "Completed."
