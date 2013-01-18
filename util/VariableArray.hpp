// VariableArray.hpp - implementation of the class used to store various
//- arrays of input information not directly used by cubit.  Normally this is a set of doubles
//- in pairs or triplets used to define a function for loadint (force and time) temperature 
//- dependence (ksi, strain%. temperature) etc.  Can also be a set of integers if needed, but
//- for now is not extensible to mixed types

// Written by Paul Wolfenbarger
// Jan 13th, 2010

#ifndef VARIABLE_ARRAY_HPP
#define VARIABLE_ARRAY_HPP
#include <vector>
#include <limits>
#include <cstddef>


template <class T>
class VariableArray : public std::vector <std::vector <T> >
{
   public :
      VariableArray (size_t rowWidth) : _rowWidth(rowWidth),
                                        mNan(std::numeric_limits<T>::quiet_NaN()) {}
      ~VariableArray() {}

      T& getItem(size_t rowNumber, size_t columnNumber) 
      {
         if (rowNumber < this->size() && columnNumber < _rowWidth)
            return (*this)[rowNumber][columnNumber] ;
         else 
            return mNan ;
      }

      void addRow(const std::vector <T>& newRow) 
      {
         this->push_back(newRow) ;
      }
      
      std::pair <int, int> getSize() // row, column
      {
         return std::make_pair(this->size(), _rowWidth) ;
      }
      
   private:
      VariableArray() ; // not implemented
      VariableArray(VariableArray const&) ; // not implemented
//       VariableArray& operator=(VariableArray const&) ; // the default should work well

      size_t _rowWidth ;
      T mNan ;
} ;

   
#define Def_Amplitude(Name) VariableArray <double> Name(2) ; 

#endif //VARIABLE_ARRAY_HPP
