#!/bin/sh


LIST=./data_5glaciers.csv
NDSI_LIST="0.45"
#NDSI_LIST="0.45 0.50 0.55"

LogDir='./log_C/'
#LogDir='./log_G/'

while read line; do
  Name=`echo ${line} | cut -d ',' -f 1`
  Lon=`echo ${line} | cut -d ',' -f 2`
  Lat=`echo ${line} | cut -d ',' -f 3`


  if [ ${Name} = 'Name' ]; then
    echo "Start"

  else
    echo "==================================================="
    echo "Glacier: "${Name}",  Lon&Lat: "${Lon}, ${Lat}
    for NDSI in ${NDSI_LIST}; do
      echo "---------------------------------------------------"
      echo "NDSI: "${NDSI}
      LOG_S2=${LogDir}${Name}_${NDSI}_S2.log
      LOG_LS8=${LogDir}${Name}_${NDSI}_LS8.log
      LOG_LS7=${LogDir}${Name}_${NDSI}_LS7.log
      LOG_LS5=${LogDir}${Name}_${NDSI}_LS5.log
      python3 Aspect_S2.py  ${Name} ${Lon} ${Lat} ${NDSI} > ${LOG_S2} &
      python3 Aspect_LS8.py ${Name} ${Lon} ${Lat} ${NDSI} > ${LOG_LS8} &
      python3 Aspect_LS7.py ${Name} ${Lon} ${Lat} ${NDSI} > ${LOG_LS7} &
      python3 Aspect_LS5.py ${Name} ${Lon} ${Lat} ${NDSI} > ${LOG_LS5} &

      NUM=`ps -U $USER | grep python3 | wc -l | awk '{print $1}'`
      while [ $NUM -gt 10 ] ; do
        sleep 30
        NUM=`ps -U $USER | grep python3 | wc -l | awk '{print $1}'`
      done

    done
  fi
done < ${LIST}

#================================================
#Finalize
NUM=`ps -U $USER | grep python3 | wc -l | awk '{print $1}'`
while [ $NUM -gt 0 ] ; do
	sleep 30
	NUM=`ps -U $USER | grep python3 | wc -l | awk '{print $1}'`
done

echo "Finished"
