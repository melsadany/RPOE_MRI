#!/bin/bash
# script will ask for XNAT credentials and project name
# then, it will download an XML with project details/participant metadata and iterate over it to download


if [ $# != 2 ] ; then
  echo "Usage: `basename $0` {XNAT_PROJ} {P_DIR}"
  exit 0;
fi


# extract input info
XNAT_PROJ=$1
P_DIR=$2

: <<'trialdata'
XNAT_PROJ="JM_GNC"
P_DIR=/Dedicated/jmichaelson-sdata-new/private-data/RPOE_MR/sourcedata/downloaded
trialdata

# get account credentials for accessing XNAT
echo -n "XNAT username: "
read username
echo -n "Password: "
read -s password

# set up connection to XNAT
exp_link="https://rpacs.iibi.uiowa.edu/xnat/data/projects/${XNAT_PROJ}/experiments?format=csv"

curl -X GET \
    -u "${username}:${password}" ${exp_link} \
    --show-error > $P_DIR/${XNAT_PROJ}.csv

# ask if you want to download a specific participant or loop
echo -n "specific participant? (Y/N)"
read sp_participant

if [ sp_participant == "N" ]; then
    while IFS="," read accessOrder id project date xsiType label insertDate URI; do
	[ "$id" == ID ] && continue;  # skips the header
	if [ ! -f ${P_DIR}/${id}.zip ]; then
            img="https://rpacs.iibi.uiowa.edu/xnat"$URI"/scans/ALL/resources/DICOM/files?format=zip"
            curl -X GET \
                -u "${username}:${password}" ${exp_link} \
                ${img} --show-error >  $P_DIR/${id}.zip
	else
            echo "File has already been downloaded for ${id}"
            continue
	fi
    done < $P_DIR/${XNAT_PROJ}.csv
else 
    echo -n "enter the participant(s) id you want (same id format from the downloaded csv id column): "
    read participant_id
    
    if [ ! -f ${P_DIR}/${participant_id}.zip ]; then
        img="https://rpacs.iibi.uiowa.edu/xnat/data/experiments/"$participant_id"/scans/ALL/resources/DICOM/files?format=zip"
        curl -X GET \
            -u "${username}:${password}" ${exp_link} \
            ${img} --show-error >  $P_DIR/${participant_id}.zip
    else
        echo "File has already been downloaded for ${participant_id}"
        continue
    fi
fi




