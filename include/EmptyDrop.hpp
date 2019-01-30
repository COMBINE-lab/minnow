#ifndef EMPTYDROP_HPP
#define EMPTYDROP_HPP

// Adopted from R version written by Lun et al. R source code 
// https://rdrr.io/github/MarioniLab/DropletUtils/src/R/emptyDrops.R
// with some changes to suit the purpose, the basic algorithm right now 
// is exactly same as shown in the code. GoodTuring test is used from the
// C++ implementation David Elworthy and devised by Sampson and Gale, 
// should be ported with this  copy of the source code. 

#include "MatrixParser.hpp"
#include <vector>

template <typename T>
void calculateEmptyDrops(DataMatrix<T>& dataMatrixObj, 
                         std::vector<bool>& emptyDropLabels) ;



#endif 