set -e

module load nco
declare -a vars=("CLDTOT" "FLUT" "FSDS" "FSDSC" \
                 "PS" "T" "TREFHT" "TREFHTMX" \  
                 "TREFHTMN" "U" "V" "Z3")

#Splits single variable timeseries into one file per year 
for i in "${vars[@]}"; do

   echo "$i"
   files=(*.$i.*)

   for file in "${files[@]}"; do

      file_name=${file%%.nc}

      year_start=${file_name: -17 }
      year_start=${year_start:0:4}
      year_end=${file_name: -8}
      year_end=${year_end:0:4}
      nyear=$(($year_end - $year_start))

      for year in $( eval echo {0..$nyear} ); do

         stride_start=$(($year * 365))
         stride_end=$(($stride_start + 364))
         year_start=${file_name: -17 }
         year_start=${year_start:0:4}
         year_current=$(($year_start + $year))

         file_out="${file::-20}"
         file_out="${file_out}${year_current}0101-${year_current}1231.nc"

         echo $file_out
         ncks -O -d time,$stride_start,$stride_end $file $file_out
      done
   done
done 
