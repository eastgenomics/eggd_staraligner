#!/bin/bash

##license related functions, different for applets
sentieon_install_dir()
{
	echo /usr/local/sentieon-genomics-*|awk '{print $1}'
}

download_license_token()
{
	#download the record contained in the $license_project
	for i in $(seq 1 5); do
		#use record instead
		RECORD_NAME=$(dx find data --class record --path $license_project --name Sentieon_License --brief|head -1)
		dx describe $RECORD_NAME --json | jq .details >~/Sentieon_License && break || sleep 5;
	done
	#cat to different file and then move, to prevent unlikely race condition
	cat ~/Sentieon_License | base64 > ~/Sentieon_License_encrypt_tmp
	mv ~/Sentieon_License_encrypt_tmp ~/Sentieon_License_encrypt
}

find_license_project()
{
	#need to select the right project incase the user belongs to multiple orgs
	billTo=$(dx describe $DX_JOB_ID --json |jq .billTo)
	launchedBy=$(dx describe $DX_JOB_ID --json|jq .launchedBy)
	#first look for the project shared with the billTo
	find_license_project_look_only $billTo
	#if not found look for all projects that the user has access to, and find the one shared with the billTo;
	#This is necessary to support the use case where the user is part of a collaborator org but not the billTo one
	if [ "$license_project" == "" ]; then
		echo "Could not find license folder containing a valid license shared with $billTo, trying to find collaborator license folders" 1>&2
		for project in $(dx api system findProjects '{"billTo":"user-sentieon_license","tags":"Sentieon_License","level":"VIEW"}'|jq .results[].id|sed 's|"||g'); do
			if [ "$project" != "null" ]; then
				for auth_user in $(dx api system findProjectMembers '{"project":"'$project'"}'|jq .results[].id|sed 's|"||g'); do
					if [ "$auth_user" == "$billTo" ]; then
						check_license_valid $project
					fi
				done
			fi
		done
	fi
	#if still not found, look for license for user
	if [ "$license_project" == "" ]; then
		echo "Could not find license folder containing a valid license shared with $billTo or a collaborator project, trying to find license folder shared with $launchedBy" 1>&2
		find_license_project_look_only $launchedBy
	fi 
	#if found expired license, set it up so that correct error message gets produced
	if [ "$license_project" == "" ] && [ "$license_project_expired" != "" ]; then
		echo "Project found but it contains an expired license"
		license_project=$license_project_expired
	fi
}

find_license_project_look_only()
{
	#if there are multiple projects, find the one that is still valid, and prioritize non-EVALs
	for project in $(dx api system findProjects '{"billTo":"user-sentieon_license","tags":"Sentieon_License","level":"VIEW","sharedWith":'$1'}'|jq .results[].id|sed 's|"||g'); do
		if [ "$project" != "null" ]; then
			check_license_valid $project
		fi
	done
}

check_license_valid()
{
	project=$1
	RECORD_NAME=$(dx find data --class record --path $project --name Sentieon_License --brief|head -1)
	RECORD=$(dx describe $RECORD_NAME --json)
	if [ $(date -d $(echo $RECORD |jq .details.license_expiration|sed 's|"||g') +"%Y%m%d") -ge $(date +"%Y%m%d") ]; then
		#replace license_project if this is purchased, to give priority to purchased
		if [ "$license_project" != "" ]; then
			if [ "$(echo $RECORD |jq .details.key|sed 's|"||g'|cut -c1-5)" != "EVAL0" ]; then
				license_project=$project
			fi
		else
			license_project=$project
		fi
	else
		license_project_expired=$project
	fi
}

error_report()
{
	dx-jobutil-report-error "$1"
	#also echo the error message to stderr, as DNAnexus does not do this
	echo "$1" 1>&2
}

