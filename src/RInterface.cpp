#include <Rcpp.h>
#include "AffineAlignment.h"
using namespace Rcpp;



Alignment Align(const std::string& target, const std::string& query, bool useAffine) {
  if (useAffine) {
    return AlignAffine(target, query);
  } else {
    return AlignSimple(target, query);
  }
}



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
List AlignRefAndRead(std::string ref, List read, bool useAffine = false) {
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

  auto aln1 = Align(ref, query, useAffine);
  auto query_rc = ReverseComplement(query);
auto aln2 = Align(ref, query_rc, useAffine);
  Alignment aln;
  if (aln1.Score >= aln2.Score) {
    aln = aln1;
    aln.QueryId = id + "/F";
  } else {
    aln = aln2;
    aln.QueryId = id + "/R";
  }
  return List::create(_[IDENTIFIERS::READ] = aln.Query,
                      _[IDENTIFIERS::REF] = aln.Ref,
                      _[IDENTIFIERS::SCORE] = aln.Score,
                      _[IDENTIFIERS::ID] = aln.QueryId);
}



// Count the number of non-gap positions in a string
int countBasePairs(const std::string& seq) {
  return seq.length() - std::count(seq.cbegin(), seq.cend(), '-');
}

// Figure out the maximum number of gaps needed in front of each reference
// basepair position.  Takes an aligned ref seq and a vector of the current
// maxes, and increases if necessary.
void addMaxGaps(const std::string& seq, std::vector<int>& sizes){
  int gapSize = 0;
  int curGapPos = 0;
  for(auto bp = seq.cbegin(); bp < seq.cend(); bp++) {
    if(*bp == '-') {
      gapSize++;
    } else {
      auto cval = sizes[curGapPos];
      if(cval < gapSize) {
        sizes[curGapPos] = gapSize;
      }
      gapSize = 0;
      curGapPos++;
    }
  }
  if(curGapPos != (sizes.size() - 1)) {
    Rcout << sizes.size() << "," << curGapPos << "\n";
    Rcpp::stop("Bug in gap finding logic.");
  }
  // Update the last pos?
  auto cval = sizes[curGapPos];
  if(cval < gapSize) {
    sizes[curGapPos] = gapSize;
  }
}

std::string gapifyReference(const std::string& ref, const std::vector<int>& gapSizes, const int fullSize) {
  std::string res;
  res.reserve(fullSize);
  auto gapPos = 0;
  int gapSize;
  for(auto it = ref.cbegin(); it < ref.cend(); it++)
  {
    auto bp = *it;
    if (bp != '-') {
      gapSize = gapSizes.at(gapPos);
      for(;gapSize > 0; gapSize--) {
        res.push_back('-');
      }
      gapPos++;
      res.push_back(bp);
    }
  }
  gapSize = gapSizes.at(gapPos);
  for(;gapSize > 0; gapSize--) {
    res.push_back('-');
  }
  return res;
}


std::string gapifySequence(const std::string& ref, const std::string& read, const std::vector<int>& gapSizes, const int fullSize) {

  std::string res;
  res.reserve(fullSize);
  // Loop through adding gaps
  int curRefBase = 0;  // Position in the ungapped reference
  int lastRefPos = -1;
  for (int curRefPos = 0; curRefPos < ref.size(); curRefPos++) {
    if (ref[curRefPos] != '-') {
      int gaps = curRefPos - (lastRefPos + 1);
      int neededGaps = gapSizes.at(curRefBase) - gaps;
      for(int i=0; i < neededGaps; i++) {
        res.push_back('-');
      }
      for(int i= (lastRefPos +1); i <= curRefPos; i++) {
        res.push_back(read.at(i));
      }
      lastRefPos = curRefPos;
      curRefBase++;
    }
  }

  // now handle the end of sequence
  int endGaps = gapSizes[gapSizes.size() - 1];
  int gapsPresent = ref.size() - (lastRefPos + 1);
  int gapsNeeded = endGaps - gapsPresent;
  for(int i=0; i< gapsNeeded; i++) {
    res.push_back('-');
  }
  for(int i=(lastRefPos+1); i < ref.size(); i++){
    res.push_back(read.at(i));
  }
  return res;
}

//' Takes a list of pairwise alignments and generates a data frame with the MSA
//' sequence. The reference sequence is determined from the first alignment in
//' the list. All reads are expanded by adding enough inserts to account for all
//' of them.
//'
//' @param ref A string that is the reference window we are aligning to.
//' @param read A list with elements name/read.
//'
//' @return Returns a list with the aligned read, ref and score
//' @export
// [[Rcpp::export]]
List AlnsToDataFrame(List alns) {

  auto ref = as<std::string>(as<List>(alns[0])[IDENTIFIERS::REF]);
  auto ref_size = countBasePairs(ref);

  // Now figure out how many "gap" positions we need in front of each reference base
  std::vector<int> gapSizes(ref_size + 1, 0); // Add before and after
  for(int i=0; i < alns.size(); i++) {
    auto cref = as<std::string>(as<List>(alns[i])[IDENTIFIERS::REF]);
    addMaxGaps(cref, gapSizes);
  }

  // Now let's modify each read so that we expand the gaps to account for each position.
  auto fullsize = ref_size + std::accumulate(gapSizes.begin(), gapSizes.end(), 0);
  auto N = alns.length() + 1;
  std::vector<std::string> seqs;
  seqs.reserve(N);
  NumericVector scores(N);
  CharacterVector ids(N);
  ids[0] = "Reference";
  scores[0] = 0;
  // Add in the reference
  seqs.push_back(gapifyReference(ref, gapSizes, fullsize));
  for(int i=0; i < alns.size(); i++)
  {
    auto elm = as<List>(alns[i]);
    auto cread = as<std::string>(elm[IDENTIFIERS::READ]);
    auto cref = as<std::string>(elm[IDENTIFIERS::REF]);
    ids[i+1] = as<std::string>(elm[IDENTIFIERS::ID]);
    scores[i+1] = elm[IDENTIFIERS::SCORE];
    seqs.push_back(gapifySequence(cref, cread, gapSizes, fullsize));
  }
  return DataFrame::create(
    _[IDENTIFIERS::ID] = ids,
    _[IDENTIFIERS::SCORE] = scores,
    _["seq"] = seqs
  );

  return List::create(_["read"] =2);
}
