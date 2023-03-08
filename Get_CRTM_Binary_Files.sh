foldername="fix_REL-3.0.0_20230303"

filename="${foldername}.tgz"
if test -f "$filename"; then
    if [ -d "fix/" ]; then #fix directory exists
        echo "fix/ already exists, doing nothing."
    else
        #untar the file and move directory to fix
				tar -zxvf $filename
				mv $foldername fix
				echo "fix/ directory created from existing $filename file."
    fi 
else
    #download, untar, move
		echo "downloading $filename, please wait about 5 minutes (3.3 GB tar file)"
    wget  ftp://ftp.ssec.wisc.edu/pub/s4/CRTM/$filename # CRTM binary files, add "-q" to suppress output. 
    tar -zxvf $filename
		mv $foldername fix
		echo "fix/ directory created from downloaded $filename."
fi
echo "Completed."
