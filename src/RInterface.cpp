#include <Rcpp.h>
#include "Alignment.h"

using namespace Rcpp;

//' Align a read and a reference window.  Both the forward and reverse direction
//' of the read will be  aligned, and only the highest scoring alignment will be
//' returned.
//'
//' @param ref A string that is the reference window we are aligning to.
//' @param read A list with elements name/read.
//'
//' @return Returns a list with the aligned read, ref and score
//' @export
// [[Rcpp::export]]
List AlignRefAndRead(std::string ref, List read) {
  // Verify read input
  if (read.size() !=2) {
    stop("Read is not a list of size 2");
  }
  CharacterVector names = read.names();
  if (names.length() != 2 || names[0] != "name" || names[1] != "read") {
    stop("Read does not have only two elements named name/read.");
  }
  auto id = as<std::string>(read["name"]);
  std::string query = as<std::string>(read["read"]);

  auto aln1 = Align(ref, query);
  auto query_rc = ReverseComplement(query);
  auto aln2 = Align(ref, query_rc);
  Alignment aln;
  if (aln1.Score >= aln2.Score) {
    aln = aln1;
    aln.QueryId = id + "/F";
  } else {
    aln = aln2;
    aln.QueryId = id + "/R";
  }
  return List::create(_["read"] = aln.Query,
                      _["ref"] = aln.Ref,
                      _["score"] = aln.Score,
                      _["id"] = aln.QueryId);
}

