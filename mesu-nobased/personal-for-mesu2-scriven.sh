#!/bin/bash
#set -x
#$position=alpha
echo working directory is $PWD 
# StandardD=./template
StandardD=../
echo standard is $StandardD

#input variables- 1
#tpvalue=(0) #Ca number
#maxlevelvalue=(64 128 256 512 1024 2048) #MESHN
#Rbubvalue=(29.8384) #Ja number
#input variables- 2
#tpvalue=(0) #Ca number
#maxlevelvalue=(256 512 1024 2048) #MESHN
#Rbubvalue=(29.8384) #Ja number
tpvalue=(15) #Ca number
#maxlevelvalue=(256 512 1024 2048) #MESHN
maxlevelvalue=(11) #(7,8,9) #MESHN:maxl-level of grid
Rbubvalue=(30) #(30) #Ja number 

# case_number='chunheng1-modify-' 
# case_number='213-basedon212-xiangbin6-myiforce-' 
# case_number='220-basedon29-xiangbin6-sigma-in-project-false-iforce2-' 
# case_number='222-basedon219-xiangbin6-sigma-in-project-false-iforcetrue-' 
# case_number='223-basedon220-xiangbin6-sigma-in-project-false-iforcetrue-cs0.1-' 
# case_number='223-basedon220-xiangbin6-sigma-in-project-false-iforcetrue-' 
# case_number='224-basedon220-xiangbin6-sigma-in-project-false-iforcetrue-jumpuv0-'
# case_number='225-new-start-multisinks-after-pre181-basedon225-5-xiangbin6-addfrom-224-' 
# case_number='225-new-start-multisinks-after-pre181-basedon225-5-xiangbin6-addfrom-224-abovevalue-225-' 
# case_number='225-new-start-multisinks-after-pre181-basedon225-5-xiangbin6-addfrom-224-abovevalue-225-Tgrad-' 
# case_number='226-new-start-multisinks-after-pre181-basedon225-5-xiangbin6-addfrom-224-abovevalue-205-Tgrad-notgood-' 
# case_number='227-new-start-multisinks-after-pre181-basedon226-sf10-' 
# case_number='231-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-'
# case_number='232-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-CFL0.1-'
# case_number='233-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-moothTgrad-' 
# case_number='234-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-smoothTgrad-noangle-' 
# case_number='235-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-smoothTgrad-noangle-tolerance-' 
# case_number='237-new-start-multisinks-after-pre181-basedon226-sf10-nodeletetriple-smoothTgrad-noangle-tolerance2-' 
case_number='244-new-start-multisinks-after-pre181-basedon237-sf10-fluxrestriction-CFL0.05-noremovedroplet-' 
for tp in "${tpvalue[@]}"; 
   do
   for rbubv in "${Rbubvalue[@]}";
    do   
      for maxlval in "${maxlevelvalue[@]}";
        do
        #title=$case_number"Ja"$rbubv"Ca"$tp"MESH"$maxlval"cpu256-ptolerance"
	title=$case_number"Ja"$rbubv"angle"$tp"MESH"$maxlval
        DestD=./$title
        rm -rf $DestD
	mkdir $DestD
	find $StandardD -maxdepth 1 -type f -exec cp {} $DestD \;
	cp -rf $StandardD"/superheated-1.0-20-5000" $DestD
       # cp -rf $StandardD $DestD
        cd $DestD
        echo now enter $PWD
        StandardIN=template01-template-2-2-restart-2-noT.c
              #template01-template-2.c
        DestIN=$title.c
       # StandardMake=template.make
        #DestMake=Makefile
        DestEXE=$title
        
        
         
      # echo $DestD $DestIN $DestMake $DestEXE
            #rm -fr input out stats testinput-* *.tmp *.root
            #let nprocs=$np_y;
            
    #sed s/CANUMBER/$tp/g $StandardIN | sed s/MESHN/$maxlval/ | sed s/JACOBNUMBER/$rbubv/ > template01.c
    sed s/CONTACTANGLE/$tp/g $StandardIN | sed s/MESHN/$maxlval/g | sed s/JACOBNUMBER/$rbubv/ > template01.c
           # sed s/CFILE/$DestIN/g $StandardMake | sed s/EXEFILE/$DestEXE/ > $DestMake
       # qcc -O2 -Wall $DestIN -o $DestEXE -lm
       sed s/MESUPROJECTNAME/$title/g ../run2-template.sh > run2.sh
       sed s/MESUPROJECTNAME/$title/g ../run2-sub-template.sh > run2-sub.sh
       echo 'here'
       #qcc -source -D_MPI=1  template01.c #-L$HOME/gl/ -lglutils -lfb_osmesa -lOSMesa -lGLU -lm
       qcc -source -grid=quadtree -D_MPI=1  template01.c
           # cd DestD
        #    make $DestMake
            sleep 3s
            #commandline=./$DestEXE
	    PASS='axSr7E-P'
	   # commandline1='sshpass -p '$PASS' ssh chenx@mesu.dsi.upmc.fr'
        
	   # commandline2='cd ~/basilisk-work ;'
	   # aimdirectory='/scratchbeta/chenx/myprogram_scratch_space3'

#NOTE!!!!!!!!!!!!!!! aimdirectory should be added in mesu, before hand
            # aimdirectory='/scratchbeta/chenx/basilisk-boiling-uniform'
       aimdirectory='/scratchbeta/chenx/axi-20231003'
	   commandline2="cd $aimdirectory ; "
	    commandline3="rm -rf $title ; mkdir $title ; "
	    commandline4="sshpass -p $PASS scp *.sh template01.c _template01.c *.h chenx@mesu.dsi.upmc.fr:~/../..$aimdirectory/$title ; sleep 5s;"
        commandline5="cd $aimdirectory/$title; chmod -R 777 run2.sh; chmod -R 777 run2-sub.sh; "
	    commandline6='./run2.sh ;'
            #xterm -e "$commandline1 && sleep 20s && logout ;" &
            #gnome-terminal "$commandline1 && sleep 10s && logout" &

USERNAME=chenx
HOSTNAME="mesu.dsi.upmc.fr"
SCRIPT1=$commandline2$commandline3
SCRIPT2=$commandline4
SCRIPT3=$commandline5$commandline6
echo $SCRIPT1
    sshpass -p $PASS ssh -l $USERNAME $HOSTNAME "$SCRIPT1"
    sleep 1s
echo $SCRIPT2
    xterm -e "$SCRIPT2" 
    sleep 20s
echo $SCRIPT3
    sshpass -p $PASS ssh -l $USERNAME $HOSTNAME "$SCRIPT3"
            sleep 15s
            rm -rf '_'$StandardIN
            #ln -s  $DestIN            
        cd ..
            echo now back to $PWD
        
    done
done

done



