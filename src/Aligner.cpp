#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

#include "Alignment.h"
using namespace Rcpp;

// Code is just a modified copy of the CC2 implementation from:
// https://github.com/PacificBiosciences/ConsensusCore2/blob/master/src/align/PairwiseAlignment.cpp


// SW parameters
const int MATCH    =  0;
const int MISMATCH = -2;
const int DELETION = -1;
const int INSERTION = -1;

inline int Max3(int a, int b, int c) { return std::max((a), std::max((b), (c))); }
inline int ArgMax3(int a, int b, int c)
{
  if (a >= b && a >= c)
    return 0;
  else if (b >= c)
    return 1;
  else
    return 2;
}

std::string Reverse(const std::string& input)
{
  std::string output;
  output.reserve(input.length());
  for (std::string::const_reverse_iterator it = input.crbegin(); it != input.crend(); ++it)
    output.push_back(*it);
  return output;
}


Alignment Align(const std::string& target, const std::string& query)
{

  int I = query.length();
  int J = target.length();
  int Score[I + 1][J + 1];

  // Fill in score matrix
  Score[0][0] = 0;
  for (int i = 1; i <= I; i++) {
    Score[i][0] = i * INSERTION;
  }
  for (int j = 1; j <= J; j++) {
    Score[0][j] = j * DELETION;
  }
  for (int i = 1; i <= I; i++) {
    for (int j = 1; j <= J; j++) {
      bool isMatch = (query[i - 1] == target[j - 1]);
      Score[i][j] = Max3(Score[i - 1][j - 1] + (isMatch ? MATCH : MISMATCH),
            Score[i - 1][j] + INSERTION, Score[i][j - 1] + DELETION);
    }
  }
  int score = Score[I][J];

  // Traceback, build up reversed aligned query, aligned target
  std::string raQuery, raTarget;
  int i = I, j = J;
  while (i > 0 || j > 0) {
    int move;
    if (i == 0) {
      move = 2;  // only deletion is possible
    } else if (j == 0) {
      move = 1;  // only insertion is possible
    } else {
      bool isMatch = (query[i - 1] == target[j - 1]);
      move = ArgMax3(Score[i - 1][j - 1] + (isMatch ? MATCH : MISMATCH),
                     Score[i - 1][j] + INSERTION, Score[i][j - 1] + DELETION);
    }
    // Incorporate:
    if (move == 0) {
      i--;
      j--;
      raQuery.push_back(query[i]);
      raTarget.push_back(target[j]);
    }
    // Insert:
    else if (move == 1) {
      i--;
      raQuery.push_back(query[i]);
      raTarget.push_back('-');
    }
    // Delete:
    else if (move == 2) {
      j--;
      raQuery.push_back('-');
      raTarget.push_back(target[j]);
    }
  }

  return Alignment("", Reverse(raQuery), Reverse(raTarget), score);
}

