In order to use Sentieon at DNAnexus you will need a license.
The code enclosed in this folder will help you use the license
you use to run the official Sentieon DNANexus apps in other
DNAnexus environments such as applets.

The folder contains 4 files:
* README.txt: this file.
* license_setup.sh: the main script that you will need to run in
  your applet code, which gathers the license info, starts the
  process that takes care of updating the license info, and sets
  up all the necessary settings for our tools to run.
* helper_funcs.sh: a library script that contains the necessary
  functions to find where the license data is located and download
  it to the worker at DNAnexus.
* license_auth.sh: a script that will take care of updating the
  license info regularly and keep it up to date.

You need to put the all 3 files (except the README.txt file) in
the home folder of your applet (/home/dnanexus) and call the
license_setup.sh script at the beginning of the applet run.
In addition, you need to make sure that the applet has external network
access, and VIEW access to your projects (to be able to find the license
information), so you will need to add the following access rules to your
dxjson.app:
 "access": { "network": ["*"], "allProjects": "VIEW" }
