#!/bin/bash
#
# Download the control FIT files dataset
#

nxsa_AI_home="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd | sed 's,/*[^/]\+/*$,,')"
dataset="${nxsa_AI_home}/dataset"
obsid_file="${nxsa_AI_home}/utils/obsids_list_test.txt"

mkdir -p $dataset
mkdir -p ${nxsa_AI_home}/tmp
cd $dataset

# Cleaning bad downloaded obs
#rm -Rf $(find . -name "data-action-aio*" -exec echo {} \; -print0 | awk -F/ '{ print $2 }')

for obs in $(cat $obsid_file | awk '{ print $1 }'); do
    if [ ! -d "${dataset}/${obs}" ]; then
        mkdir -p ${dataset}/${obs}/pps
        cd ${dataset}/${obs}/pps
        echo "***  Downloading data for obsid: ${obs} ..."

        # Events FIT file for PN camera
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=PIEVLI&level=PPS&instname=PN&extension=FTZ" > ${nxsa_AI_home}/tmp/${obs}_events.log 2>&1 &

        # Regions for all cameras
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=REGION&level=PPS&extension=ASC" > ${nxsa_AI_home}/tmp/${obs}_regions.log 2>&1 &

        # Image FIT file for PN camera
        curl -L -O -J "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obs}&name=IMAGE_8&level=PPS&instname=PN&extension=FTZ" > ${nxsa_AI_home}/tmp/${obs}_image.log 2>&1 &

        wait
    else
        echo "*** Folder ${obs} exists ..."
    fi

done
