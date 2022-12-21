#!/bin/bash
#
license_project=$1

source ~/helper_funcs.sh

#refresh the token file
while sleep 300; do
	download_license_token
done
