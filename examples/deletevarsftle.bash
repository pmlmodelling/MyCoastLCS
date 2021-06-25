for ff in $( ls *nc )
do
	fileout=${ff:(-3)} # selects the last 3 characters
        fileout=${ff:1:15} # selects the first 15 character
        fileout=${ff#.nc} # hopefully removes '.nc' from the filename starting at the front
        fileout=${ff%.nc} # this should start from the back 	
	echo ${fileout}_stripped.nc
	ncks -x -v x,y,z $ff ${fileout}_stripped.nc
       echo "Done removing the big variables" 
       echo "Now I should remove the times that are identical..."

done
