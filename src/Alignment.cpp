#include "Alignment.h"
Alignment::Alignment(std::string id, std::string q, std::string r, int score) :
  QueryId(id),
  Query(q),
  Ref(r),
  Score(score)
{}

std::string ReverseComplement(const std::string& seq) {
  std::string output;
  output.reserve(seq.length());
  for (auto it = seq.crbegin(); it != seq.crend(); ++it) {
    auto bp = *it;
    switch(bp) {
      case 'A':
        output.push_back('T');
        break;
      case 'C':
        output.push_back('G');
        break;
      case 'G':
        output.push_back('C');
        break;
      case 'T':
        output.push_back('A');
        break;
      case '-':
        output.push_back('-');
        break;
      default:
        Rcpp::stop("Tried to reverse complement sequence with non ATGC-");
    }
  }
  return output;
}
