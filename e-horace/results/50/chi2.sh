pseudofile=oldalpha/tev_z_mu_minv_oal_norec_200.dat
mz=90687

for ((i=1; i<=21; ++i))
do

    teorfile=born/$i/tev_z_mu_minv_born_norec_200.dat
    chi=0
    
    for ((line=61; line<=131; ++line))
    do
    
        
        
        tbin=$(sed -n $line'p' $teorfile | awk '{print $2}')
        if [[ $tbin =~ [Ee] ]];then
        
           base1=$(echo $tbin | cut -d 'E' -f1)
           exp1=$(($(echo $tbin | cut -d 'E' -f2)*1))
           tbin=$(bc -l <<< "$base1*(10^$exp1)")
        fi
        
               
        pbin=$(sed -n $line'p' $pseudofile | awk '{print $2}')
        if [[ $pbin =~ [Ee] ]];then
        
           base2=$(echo $pbin | cut -d 'E' -f1)
           exp2=$(($(echo $pbin | cut -d 'E' -f2)*1))
           pbin=$(bc -l <<< "$base2*(10^$exp2)")
        fi
        
              
        err=$(sed -n $line'p' $pseudofile | awk '{print $3}')
        base=$(echo $err | cut -d 'E' -f1)
        exp=$(($(echo $err | cut -d 'E' -f2)*1))
        err=$(bc -l <<< "$base*(10^$exp)")

        chi=$(echo "scale=17; $chi+($tbin-$pbin)*($tbin-$pbin)/($err*$err)" | bc) 

    
        

    done 
    
   
    echo $mz $chi  >> chioldalpha.dat
    
    mz=$((mz+50))

done
 

exit
