///
/// \file     coo_reader.hpp
/// \author   zhouxingbin
/// \group    CDCS-HPC
/// \date     2022-1-19
/// \brief    COO矩阵读取。
///
#ifndef __CHIPSUM_COO_READER_HPP__
#define __CHIPSUM_COO_READER_HPP__

#include <fstream>
#include <sstream>



namespace ChipSum {
namespace Common {

template <typename Stype, typename Vtype>
void coo_reader(Stype &nrows, Stype &ncols, Stype &nnz,
                Stype* &row_map, Stype* &col_map, Vtype* &values, const char *fileName) {
  //read file
  std::ifstream file(fileName, std::ifstream::in);
  if (!file.is_open()) {
    throw std::runtime_error ("File cannot be opened\n");
  }
  
  std::string fline = "";
  getline(file, fline);
  
  //check first line
  if (fline.size() < 2 || fline[0] != '%' || fline[1] != '%'){
    throw std::runtime_error ("Invalid file. Line-1\n");
  }
  if (fline.find("matrix") == std::string::npos || fline.find("coordinate") == std::string::npos){
    throw std::runtime_error ("Invalid file. Line-1\n");
  }
  //skip comments
  while(1){
  getline(file, fline);
  if(fline[0] != '%') break;
  }
  //read coo meta data
  std::stringstream ss (fline);
  ss >> nrows >> ncols >> nnz;
  
  //set memory
  row_map = new Stype[nnz];
  col_map = new Stype[nnz];
  values = new Vtype[nnz];

  for (Stype i = 0; i < nnz; ++i){
    getline(file, fline);
    std::stringstream ss2(fline);
    ss2 >> row_map[i] >> col_map[i] >> values[i];  
  }

    
file.close();
}

}
}
#endif
