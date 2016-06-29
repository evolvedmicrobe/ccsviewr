#include <string>
#include <Rcpp.h>

// Namespace to hold strings that will refer to different list elements that
// are passed back and forth between C++ and R.
namespace IDENTIFIERS {
  const std::string ID = "id";
  const std::string READ = "read";
  const std::string REF = "ref";
  const std::string SCORE = "score";
}

// Basic alignment class, returned by the Align method.
class Alignment {
  public:
    std::string QueryId;
    std::string Query;
    std::string Ref;
    int Score;
    Alignment(std::string id, std::string q, std::string r, int Score);
    Alignment() = default;
};

std::string Reverse(const std::string& input);

Alignment AlignSimple(const std::string& target, const std::string& query);

std::string ReverseComplement(const std::string& seq);
