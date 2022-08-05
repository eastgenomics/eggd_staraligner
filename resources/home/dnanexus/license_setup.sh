#!/bin/bash

source ~/helper_funcs.sh
SENTIEON_INSTALL_DIR=$(sentieon_install_dir)
sentieon_procs=$(nproc)
# Setup license
# find corresponding project where license token information is and download it
find_license_project
if [ "$license_project" == "" ]; then
	error_report "You do not have a license to use this app. Please contact Sentieon at info@sentieon.com to request a license. The app will exit now."
fi
download_license_token

# Kick token refresh script
bash ~/license_auth.sh $license_project&

# Set license settings
job_tag=$(cat ~/Sentieon_License|jq '.job_tag')
license_server_location=$(cat ~/Sentieon_License|jq '.license_server_location'|sed 's|"||g')
if [ "$license_server_location" == "null" ] || [ "$license_server_location" == "" ]; then
	#try global location
	set +e
	license_server_location=$(curl -s https://sentieon-bundle.s3.amazonaws.com/dnanexus/DNAnexus_app|jq '.license_server_location'|sed 's|"||g')
	set -e
	if [ "$license_server_location" == "null" ] || [ "$license_server_location" == "" ]; then
		license_server_location=master.sentieon.com:9010
	fi
fi
export SENTIEON_LICENSE=$license_server_location
export SENTIEON_AUTH_MECH=dnanexus_app
export SENTIEON_AUTH_DATA=~/Sentieon_License_encrypt
export SENTIEON_JOB_TAG=$job_tag

# Error checking: Check license validity before running anything
${SENTIEON_INSTALL_DIR}/bin/sentieon licclnt ping -s $SENTIEON_LICENSE 2> >(tee lic_errlog) || error_report "There is an issue with your license. Please contact Sentieon at support@sentieon.com and report your username, org, the billing org you are using to run Sentieon, and the following error message \"$(cat lic_errlog)\". The app will exit now."

#other settings
export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so.1
export MALLOC_CONF=lg_dirty_mult:-1

