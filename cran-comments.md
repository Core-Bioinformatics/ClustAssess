# ClustAssess v1.0.0

## Test environments
- local Ubuntu 22.04, R 4.4.0
- win-builder (oldrelease, release, devel)
- mac-builder (release, development)

## R CMD check results
There were no ERRRORs or WARNINGs.

### On local Ubuntu, there were 4 NOTESs:
> checking CRAN incoming feasibility ... [5s/65s] NOTE

Maintainer: ‘Andi Munteanu <am3019@cam.ac.uk>’

New maintainer:
  Andi Munteanu <am3019@cam.ac.uk>
Old maintainer(s):
  Arash Shahsavari <as3006@cam.ac.uk>

Suggests or Enhances not in mainstream repositories:
  monocle3

Package has a VignetteBuilder field but no prebuilt vignette index.

> checking package dependencies ... NOTE

Imports includes 32 non-default packages.

Importing from so many packages makes the package vulnerable to any of
them becoming unavailable.  Move as many as possible to Suggests and
use conditionally.

> checking installed package size ... NOTE

installed size is  7.4Mb

sub-directories of 1Mb or more:
R      1.1Mb
libs   3.9Mb

> checking for future file timestamps ... NOTE
unable to verify current time

*comment* : this NOTE is not related to the functionality of the package.

For the other environments, the NOTEs were the same as the ones above.