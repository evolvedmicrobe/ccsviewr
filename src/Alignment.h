#include <string>
#include <Rcpp.h>
class Alignment {
  public:
    std::string QueryId;
    std::string Query;
    std::string Ref;
    int Score;
    Alignment(std::string id, std::string q, std::string r, int Score);
    Alignment() = default;
};

Alignment Align(const std::string& target, const std::string& query);
std::string ReverseComplement(const std::string& seq);
