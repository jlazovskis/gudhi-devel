#ifndef SIRUP_METHODS_H
#define SIRUP_METHODS_H

#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistence_matrix/RU_matrix.h>
#include <gudhi/Persistence_matrix/ru_pairing.h>
#include <gudhi/Persistence_matrix/columns/entry_types.h>

namespace Gudhi {
namespace persistence_matrix {

template <class Matrix_entry, typename Index>
static Index get_row_index_shifted(Index row_index, auto it_begin, auto it_end) {
  while (it_begin != it_end && *it_begin < row_index) {
    row_index--;
    it_begin++;
  }
  return row_index;
}

template <class Master_matrix>
class SiRUP_methods {
  using Column = typename Master_matrix::Column;
  using Index = typename Master_matrix::Index;

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;

  // void _update_barcode(Index birthPivot, Index death) {
  //   if constexpr (Master_matrix::Option_list::has_column_pairings) {
  //     typename RU_pairing<Master_matrix>::_update_barcode(birthPivot, death);
  //   }
  // }

 public:
  static void remove_maximal_cell(Master_RU_matrix& RU, std::vector<Index> removeColumns, bool is_index_updated) {
    static_assert(Master_matrix::Option_list::has_removable_columns,
                  "'remove_maximal_cell' is not implemented for the chosen options.");

    // auto op = RU.operators_;
    // auto R = RU.reducedMatrixR_;
    // auto U = RU.mirrorMatrixU_;

    // bool iIsPositive = _matrix()->reducedMatrixR_.is_zero_column(0);

    // std::cout << "--------------------\n";
    // std::cout << "**barcode size = " << R->_birth(0) << "\n";
    // std::cout << "--------------------\n";

    // Container for columns already removed. Indexation does not change.
    std::set<Index> removedColumns;

    // Steps of algorithm implementation
    _execute_sirup(RU, removeColumns, removedColumns, is_index_updated);
    _clear_rows_columns(RU, removedColumns, is_index_updated);
    _shift_persistence_indices(RU, removedColumns);

    // R.print();
    // U.print();
  }

  // Step 1: Execute SiRUP
  static void _execute_sirup(Master_RU_matrix& RU, std::vector<Index> removeColumns, std::set<Index>& removedColumns,
                             bool is_index_updated) {
    // Shorthand
    auto op = RU.operators_;
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;
    using M = Master_matrix;
    constexpr const Index nullDeath = M::template get_null_value<Index>();

    // Iterate over columns requested to remove
    for (int i = 0; i < removeColumns.size(); i++) {

      // Reindex if necessary
      Index columnIndex = removeColumns[i];
      if (is_index_updated && removedColumns.size() > 0) {
        for (auto it_rc = removedColumns.begin(); it_rc != removedColumns.end(); ++it_rc) {
          if (*it_rc <= columnIndex) {
            columnIndex += 1;
          }
        }
      }
      // std::cout << "removing at index " << columnIndex << std::endl;

      // Get row of V (from its inverse U)
      // Note this differs from "A" in paper, as this A contains diagonal element
      // const auto A = _get_inverse_row(RU, columnIndex, removedColumns);

      const Column A_const = Column(U.get_column(columnIndex), nullptr);
      Column A = Column(U.get_column(columnIndex), nullptr);
      auto itAc = std::next(A_const.begin());
      while (itAc != A_const.end()) {
        Column cur_row = Column(U.get_column(*itAc), nullptr);
        cur_row.clear(*itAc);
        cur_row *= op->get_characteristic() - M::get_element(*itAc);
        A += cur_row;
        ++itAc;
      }

      // std::cout << "* Affected simplices: ";
      // for ( auto a : A ) { std::cout << a << " "; }
      // std::cout << std::endl;
      // R.print();

      // Construct vector of earliest-highest pivot columns
      std::vector<Index> B = {columnIndex};
      auto itA = A.begin();
      itA++;
      while (itA != A.end() && !R.is_zero_column(M::get_row_index(*itA))) {
        // while (itA != A.end() && !R.is_zero_column(row_shift(M::get_row_index(*itA)))) {
        Index b = R.get_pivot(M::get_row_index(*itA));
        // Index b = R.get_pivot(row_shift(M::get_row_index(*itA)));
        // std::cout << "** col " << *itA << " has pivot " << b << "\n";
        if (b < R.get_pivot(B.back())) {
          B.push_back(M::get_row_index(*itA));
          // B.push_back(row_shift(M::get_row_index(*itA)));
        }
        itA++;
      }
      if (itA != A.end()) {
        B.push_back(M::get_row_index(*itA));
        // B.push_back(row_shift(M::get_row_index(*itA)));
      }
      // std::cout << "* Pareto front: ";
      // for ( auto b : B ) { std::cout << b << " "; }
      // std::cout << std::endl;

      // Perform column additions and barcode reindexation
      // If column has not been added to other columns, only update persistence
      if (A.size() == 1) {
        if (!R.is_zero_column(columnIndex)){
          auto& itBar = RU.indexToBar_.at(columnIndex);
          itBar->death = nullDeath;
          RU.indexToBar_.try_emplace(itBar->birth, itBar); // Ask Hannah: is this necessary?
        } else {
          RU.barcode_.erase(RU.indexToBar_[columnIndex]);
        }
      // If column has been added to other columns
      } else {
        itA = A.end();
        itA--;
        bool targetInB = false;
        bool targetIsZeroCol = false;
        while (itA != A.begin()) {
          auto itB = B.begin();
          // std::cout << "*considering col " << M::get_row_index(*itA);
          // if (R.is_zero_column(M::get_row_index(*itA))) { std:: cout << " (is zero)\n";}
          // else { std:: cout << " (is not zero)\n"; }

          // Separate case for when target is zero column
          if (R.is_zero_column(M::get_row_index(*itA))) {
            itB = std::prev(B.end());
            if (*itB == M::get_row_index(*itA)) {
              itB--;
              targetInB = true;
            }
            targetIsZeroCol = true;
          }
          // Separate case for nonzero columns
          else {
            Index k = R.get_pivot(M::get_row_index(*itA));
            while ((R.get_pivot(*itB) > k) && (R.get_pivot(*itB) != -1)) {
              itB++;
            }
            if (*itB == M::get_row_index(*itA)) {
              itB--;
              targetInB = true;
            }
          }

          // Update bars
          if (targetInB) {
            // std::cout << "*now birth " << R.get_pivot(*itB) << " has death " << M::get_row_index(*itA) << std::endl;

            // Update for addition target, when target is last element of B
            if (std::next(itB) == std::prev(B.end())) {
              // Case 1: target is zero column, so remove bar
              if (targetIsZeroCol) {
                RU.barcode_.erase(RU.indexToBar_[M::get_row_index(*itA)]);
                // RU.indexToBar_.erase(M::get_row_index(*itA));
              // Case 2: target is nonzero column, so update bar to infinite bar
              } else {
                auto& itBar = RU.indexToBar_.at(M::get_row_index(*itA));
                itBar->death = nullDeath;
                // std::cout << "** last B update: " << *itBar << std::endl;
                // RU.indexToBar_.try_emplace(itBar->birth, itBar); // Ask Hannah: is this necessary?
              }
            }

            // Update for addition source
            RU._update_barcode(R.get_pivot(*itB), M::get_row_index(*itA));
            auto& itBar = RU.indexToBar_.at(R.get_pivot(*itB));
            RU.indexToBar_[M::get_row_index(*itA)] = itBar;
            targetInB = false;
          }

          // Perform addition
          // std::cout << "*adding " << *itB << " to " << M::get_row_index(*itA) << std::endl;
          R.add_to(*itB, M::get_row_index(*itA));
          U.add_to(M::get_row_index(*itA), *itB);
          itA--;
        }
      }

      // Ensure columns are zero
      R.zero_column(columnIndex);
      U.zero_column(columnIndex);
      removedColumns.insert(columnIndex);

      // Remove reference to bar for removed column
      RU.indexToBar_.erase(columnIndex);

      // R.print();
      // U.print();
    }
  };

  // Part of Step 1: Get inverse row
  static Column _get_inverse_row(Master_RU_matrix& RU, Index r, std::set<Index>& removedColumns) {
    // Shorthand
    auto op = RU.operators_;
    auto U = RU.mirrorMatrixU_;
    using M = Master_matrix;
    using E = typename M::Element;

    // Definitions
    const Column res_const = Column(U.get_column(r), nullptr);
    Column res = Column(U.get_column(r), nullptr);
    auto size = U.get_number_of_columns();
    auto it = std::next(res_const.begin());

    // Construct inverse
    while (it != res_const.end()) {
      Column cur_row = Column(U.get_column(*it), nullptr);
      cur_row.clear(*it);
      cur_row *= op->get_characteristic() - M::get_element(*it);
      res += cur_row;
      ++it;
    }

    // Return result
    return res;
  };

  // Step 2: Clear the rows and columns
  static void _clear_rows_columns(Master_RU_matrix& RU, std::set<Index>& removedColumns, bool is_index_updated) {
    // Remove columns
    // R.erase_empty_column(columnIndex);
    // U.erase_empty_column(columnIndex);

    auto op = RU.operators_;
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;

    auto size = R.get_number_of_columns();
    int shift = 0;

    // Swap all rows in one sweep
    if (removedColumns.size() > 0) {
      auto it_rc = removedColumns.begin();
      bool at_rc_end = false;
      for (int i = *it_rc; i < size; i++) {
        if (!at_rc_end && *it_rc == i) {
          // std::cout << "*index " << i << " was dropped" << "\n";
          shift += 1;
          it_rc++;
          if (it_rc == removedColumns.end()) {
            at_rc_end = true;
          }
        } else {
          U.swap_rows(i - shift, i);
          // std::cout << "*swapping " << i - shift << " and " << i << "\n";
        }
      }
    }

    // Erase all columns and rows in one sweep
    shift = 0;
    for (Index columnIndex : removedColumns) {
      R.erase_empty_column(columnIndex - shift);
      U.erase_empty_column(columnIndex - shift);
      U.erase_empty_row(size - 1 - shift);
      ++shift;
    }

    // R.matrix_.erase(R.matrix_.begin()+columnIndex);
    // U.matrix_.erase(U.matrix_.begin()+columnIndex);
    // U.erase_empty_row(columnIndex);

    // R.Column_container.erase(columnIndex);
    // U.Column_container.erase(columnIndex);

    // R.erase_empty_column(size-1);
    // U.erase_empty_column(size-1);
    // U.erase_empty_row(size-1);
    // R.print();
    // U.print();
  };

  // Step 3: Update barcode and index dictionary
  static void _shift_persistence_indices(Master_RU_matrix& RU, std::set<Index>& removedColumns){

    auto op = RU.operators_;
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;

    auto size = R.get_number_of_columns();
    int shift = 0;

    if (removedColumns.size() > 0) {
      auto it_rc = removedColumns.begin();
      bool at_rc_end = false;
      for (int i = *it_rc; i < size; i++) {
        if (!at_rc_end && *it_rc == i) {
          shift += 1;
          it_rc++;
          if (it_rc == removedColumns.end()) {
            at_rc_end = true;
          }
        } else {

          // Update barcode_
          auto& itBar = RU.indexToBar_.at(i);
          if (i == itBar->birth){
            // std::cout << "*column " << i << " is positive simplex" << "\n";
            itBar->birth -= shift;
          } else {
            // std::cout << "*column " << i << " is negative simplex" << "\n";
            itBar->death -= shift;          
          }
  
          // Update indexToBar_
          RU.indexToBar_.erase(i);
          RU.indexToBar_.try_emplace(i-shift, itBar); // Ask Hannah: is this necessary?
        }
      }      
    }

    // Update 'Barcode barcode_;'
    // Is of type 'std::list<Bar>'
    //  using Bar = Persistence_interval<Dimension, Pos_index>;
    //  using Pos_index = typename PersistenceMatrixOptions::Index; 
    //  Persistence_interval is a stuct with dim, birth, death

    // Update 'Dictionary indexToBar_;'
    // Is of type 'std::unordered_map<Pos_index, typename Barcode::iterator>'

    // for (const auto & [ key, value ] : RU.indexToBar_) {
    //   std::cout << key << " : " << (*value).dim << " / " << (*value).birth << " / " << (*value).death << std::endl;
    // }

  };
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // SIRUP_METHODS_H
