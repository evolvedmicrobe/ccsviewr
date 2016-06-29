#include <Rcpp.h>
using namespace Rcpp;

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <string>
#include <vector>

#include "AffineAlignment.h"

class IupacAware;
class Standard;

static float MAX4(float a, float b, float c, float d)
{
  return std::max(std::max(a, b), std::max(c, d));
}




inline float MatchScore(char t, char q, float matchScore, float mismatchScore,
                                  float partialMatchScore)
{
  return (t == q ? matchScore : mismatchScore);
}


Alignment AlignAffine(const std::string& target, const std::string& query,
                                      AffineAlignmentParams params)
{
  // Implementation follows the textbook "two-state" affine gap model
  // description from Durbin et. al

  int I = query.length();
  int J = target.length();
  NumericMatrix M(I + 1, J + 1);
  NumericMatrix GAP(I + 1, J + 1);

  // Initialization
  M(0, 0) = 0;
  GAP(0, 0) = -FLT_MAX;
  for (int i = 1; i <= I; ++i) {
    M(i, 0) = -FLT_MAX;
    GAP(i, 0) = params.GapOpen + (i - 1) * params.GapExtend;
  }
  for (int j = 1; j <= J; ++j) {
    M(0, j) = -FLT_MAX;
    GAP(0, j) = params.GapOpen + (j - 1) * params.GapExtend;
  }

  // Main part of the recursion
  for (int i = 1; i <= I; ++i) {
    for (int j = 1; j <= J; ++j) {
      float matchScore = MatchScore(target[j - 1], query[i - 1], params.MatchScore,
                                       params.MismatchScore, params.PartialMatchScore);
      M(i, j) = std::max(M(i - 1, j - 1), GAP(i - 1, j - 1)) + matchScore;
      GAP(i, j) = MAX4(M(i, j - 1) + params.GapOpen, GAP(i, j - 1) + params.GapExtend,
          M(i - 1, j) + params.GapOpen, GAP(i - 1, j) + params.GapExtend);
    }
  }

  // Perform the traceback
  const int MATCH_MATRIX = 1;
  const int GAP_MATRIX = 2;

  std::string raQuery, raTarget;
  int i = I, j = J;
  int mat;
  double score;
  if ( M(I, J) >= GAP(I, J) ) {
    mat = MATCH_MATRIX;
    score = M(I, J);
  } else {
    mat = GAP_MATRIX;
    score = GAP(I, J);
  }
  int iPrev, jPrev, matPrev;
  while (i > 0 || j > 0) {
    if (mat == MATCH_MATRIX) {
      matPrev = (M(i - 1, j - 1) >= GAP(i - 1, j - 1) ? MATCH_MATRIX : GAP_MATRIX);
      iPrev = i - 1;
      jPrev = j - 1;
      raQuery.push_back(query[iPrev]);
      raTarget.push_back(target[jPrev]);
    } else {
      float s[4];
      s[0] = (j > 0 ? M(i, j - 1) + params.GapOpen : -FLT_MAX);
      s[1] = (j > 0 ? GAP(i, j - 1) + params.GapExtend : -FLT_MAX);
      s[2] = (i > 0 ? M(i - 1, j) + params.GapOpen : -FLT_MAX);
      s[3] = (i > 0 ? GAP(i - 1, j) + params.GapExtend : -FLT_MAX);
      int argMax = std::max_element(s, s + 4) - s;

      matPrev = ((argMax == 0 || argMax == 2) ? MATCH_MATRIX : GAP_MATRIX);
      if (argMax == 0 || argMax == 1) {
        iPrev = i;
        jPrev = j - 1;
        raQuery.push_back('-');
        raTarget.push_back(target[jPrev]);
      } else {
        iPrev = i - 1;
        jPrev = j;
        raQuery.push_back(query[iPrev]);
        raTarget.push_back('-');
      }
    }

    // Go to previous square
    i = iPrev;
    j = jPrev;
    mat = matPrev;
  }
  return Alignment("", Reverse(raQuery), Reverse(raTarget), static_cast<int>(score));
}


AffineAlignmentParams::AffineAlignmentParams(float matchScore, float mismatchScore, float gapOpen,
                                             float gapExtend, float partialMatchScore)
  : MatchScore(matchScore)
, MismatchScore(mismatchScore)
, GapOpen(gapOpen)
, GapExtend(gapExtend)
, PartialMatchScore(partialMatchScore)
{
}

AffineAlignmentParams DefaultAffineAlignmentParams()
{
  return AffineAlignmentParams(0, -1.0, -1.0, -0.5, 0);
}



