for i in *.csv
do
 echo "${i}"
 awk -F "\"*,\"*" '{print $3,$4}' "${i}"
done
