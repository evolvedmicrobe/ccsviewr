#' Get an alignment of a CCS read to the reference as well as all of it's subreads
#' aligned to the reference.  This function performs the following steps:
#'
#'  1 - Gets the region of the reference genome that the CCS read aligns to
#'  2 - Gets the CCS read and the subreads
#'  3 - Aligns all reads to the section of the reference genome
#'  4 - Returns a list of all these pairwise alignments with the CCS read alignment
#'  as the first element.
#'
#'
#' @param hole Which ZMW to collect data for.
#' @param ccs_name Full path of the aligned CCS BAM file.
#' @param subreads_name Full path of the subreads BAM file
#' @param fasta_name Full path of the fasta file used to generate the CCS alignment
#'
#' @return Returns a list of alignments with the CCS read first and all subreads next
#' @export
getAlignments <- function(hole, ccs_name, subreads_name, fasta_name) {

  # Load the CCS read and the reference
  ind = pbbamr::loadpbi(ccs_name)
  cur = ind[ind$hole==hole,]
  if(nrow(cur) != 1) {
    stop("Could not find hole in CCS BAM index file.")
  }

  ref = pbbamr::loadReferenceWindow(as.character(cur$ref), cur$tstart, cur$tend, fasta_name)
  ccs_seq = pbbamr::loadSubreadsAtOffsets(cur$offset, ccs_name)[[1]]
  ccs_aln = AlignRefAndRead(ref, ccs_seq)

  # Now get the subreads
  ind_sub = pbbamr::loadpbi(subreads_name)
  subs = ind_sub[ind_sub$hole==hole,]
  if(nrow(subs)==0) {
    stop("No subreads found for the hole in the BAM file.")
  }
  subreads = pbbamr::loadSubreadsAtOffsets(subs$offset, subreads_name)
  names = sapply(subreads, function(x) x$name)
  alns = lapply(subreads, function(z) AlignRefAndRead(ref, z))
  c(list(ccs_aln), alns)
}



#' Determine the reference position in a multiple sequence alignment returned
#' from AlnsToDataFrame.
#'
#' Typically, one will call plotMSA to view a variant call in an alignment, but
#' that function takes coordinates in terms of alignment positions, not
#' reference genome positions.  This function converts a reference genome
#' location to an alignment location.  It does this by taking the start of the
#' alignment on the genome, and then figures out how much further it needs to go
#' to get the respective genome position in the alignment.
#'
#' @param tstart The initial position of the alignment.
#' @param pos The
#' @param df Data frame returned by AlnsToDataFrame.  The seq in the
#' @param fasta_name Full path of the fasta file used to generate the CCS alignment
#'
#' @return Returns The genome position
#' @export
getRefPosition <- function(tstart, pos, df) {
  neededPos = pos - tstart
  ref = strsplit(as.character(df$seq[1]), "")[[1]]
  poses = which(ref%in%c("A","C","G","T"))
  poses[neededPos]
}
