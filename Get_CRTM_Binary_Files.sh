#https://bin.ssec.wisc.edu/pub/s4/CRTM/file/crtm_coefficients_3.1.0_skylab_7.0.tar.gz
foldername="3.1.0_skylab_7.0"
filename="crtm_coefficients_${foldername}.tar.gz"
echo "$filename"
break

if test -f "$filename"; then
    if [ -d "fix/" ]; then #fix directory exists
        echo "fix/ already exists, doing nothing."
    else
        #untar the file and move directory to fix
				tar -zxvf $filename
				mkdir fix
				mv crtm/$foldername/* fix/.
				rm -rf $foldername
				echo "fix/ directory created from existing $filename file."
    fi 
else
    #download, untar, move
		echo "Downloading $filename, please wait about 5 minutes (4 GB tar file)"
	  wget  https://bin.ssec.wisc.edu/pub/s4/CRTM/file/$filename # CRTM binary files, add "-q" to suppress output. 
		
		
    #untar the file and move directory to fix
    tar -zxvf $filename
    mkdir fix
    mv crtm/$foldername/* fix/.
    rm -rf crtm/$foldername
  	echo "fix/ directory created from downloaded $filename."
fi
echo "Completed."
