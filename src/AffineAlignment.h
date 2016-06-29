//
// Support for pairwise alignment with an affine gap penalty.
//

#pragma once

#include <string>
#include <vector>
#include "Alignment.h"

struct AffineAlignmentParams
{
  float MatchScore;
  float MismatchScore;
  float GapOpen;
  float GapExtend;
  float PartialMatchScore;

  AffineAlignmentParams(float matchScore, float mismatchScore, float gapOpen, float gapExtend,
                        float partialMatchScore = 0);
};

AffineAlignmentParams DefaultAffineAlignmentParams();

//
// Affine gap-penalty alignment.
//
Alignment AlignAffine(
    const std::string& target, const std::string& query,
    AffineAlignmentParams params = DefaultAffineAlignmentParams());  // NOLINT
