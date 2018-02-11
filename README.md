#ccsviewr - PDF images of CCS consensus sequences and subreads

An R package to produce PDFs of alignments of CCS reads to a reference template along with every subread that was used to produce that CCS read.  This package aims to allow for users to gain certainty in a SNP call or to investigate problems with consensus accuracy.

An example sequence alignment shown below:


![Alignment Image](https://github.com/PacificBiosciences/ccsviewr/blob/master/vignettes/AlignmentWindow.png)

## Documentation


[The online vignette](http://htmlpreview.github.io/?http://github.com/PacificBiosciences/ccsviewr/blob/master/vignettes/ccsviewr.html)

### Implementation details

R Package Dependencies

	- pbbamr
	- grid
	- Rcpp
	
The code uses grid to plot the multiple sequence alignments, and prints to a PDF with Courier font by default.  The aligner used is a simple Smith Waterman, but with a deletion/insertion couplet more likely than a mismatch in order to mirror the PacBio error modes.

DISCLAIMER
----------
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
