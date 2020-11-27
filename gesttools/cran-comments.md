## Test environments
- Windows-x86_64-w64-mingw32/x64 (64-bit)-local install, R version 4.0.3
- R-hub windows-x86_64-devel (r-devel)
- R-hub ubuntu-gcc-release (r-release)
- R-hub fedora-clang-devel (r-devel)

## R CMD check results
> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Daniel Tompsett <danieltompsettwork@gmail.com>'
  
  New submission
  
  Possibly mis-spelled words in DESCRIPTION:
    Vansteelandt (8:126, 8:160)
    Sjolander (8:112)
  Found the following (possibly) invalid DOIs:
    DOI: 10.1093/aje/kwx347
      From: DESCRIPTION
      Status: libcurl error code 56:
      	Send failure: Connection was reset
      Message: Error

> On windows-x86_64-devel (r-devel), ubuntu-gcc-release (r-release), fedora-clang-devel (r-devel)
  checking examples ... NOTE
  Examples with CPU (user + system) or elapsed time > 5s
               user system elapsed
  gestmultcat 12.23   0.19   12.47
  gestmult     9.04   0.39    9.50
  gestcat      7.98   0.19    8.22
  gest         5.95   0.17    6.15
  gest.boot    5.50   0.00    5.52

0 errors √ | 0 warnings √ | 2 notes x

An additional note. This is a new submission by the user

The possibly mis-spelt words are the authors
and are correctly spelt as given in the referenced paper. 

The dois that are given are correct, and given in the correct format 
to our knowledge. We are unsure of the reason for the libcurl error.

