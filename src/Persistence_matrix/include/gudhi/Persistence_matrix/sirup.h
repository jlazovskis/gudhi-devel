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
    _update_persistence(RU, removedColumns);

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

    // Iterate over columns requested to remove
    for (int i = 0; i < removeColumns.size(); i++) {
      Index columnIndex = removeColumns[i];
      if (is_index_updated && removedColumns.size() > 0) {
        for (auto it_rc = removedColumns.begin(); it_rc != removedColumns.end(); ++it_rc) {
          if (*it_rc <= columnIndex) {
            columnIndex += 1;
          }
        }
      }
      // std::cout << "removing at index " << columnIndex << std::endl;

      // Definitions
      Index size = R.get_number_of_columns();
      // std::cout << "* size: " << size <<"\n";

      // Get row of V (by inverting U)
      const auto A = get_inverse_row(RU, columnIndex, removedColumns);
      // std::cout << "* Affected simplices: ";
      // for ( auto a : A ) { std::cout << a << " "; }
      // std::cout << std::endl;
      // R.print();

      // construct vector of earliest-highest pivot columns
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

      // Perform column additions
      itA = A.end();
      itA--;
      bool targetInB = false;
      while (itA != A.begin()) {
        auto itB = B.begin();
        // std::cout << "*considering col " << M::get_row_index(*itA);
        // if (R.is_zero_column(M::get_row_index(*itA))) { std:: cout << " (is zero)\n";}
        // else { std:: cout << " (is not zero)\n"; }

        // Separate case for zero columns
        if (R.is_zero_column(M::get_row_index(*itA))) {
          itB = std::prev(B.end());
          if (*itB == M::get_row_index(*itA)) {
            itB--;
            targetInB = true;
          }
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

        // Update barcode
        if (targetInB) {
          std::cout << "*now birth " << R.get_pivot(*itB) << " has death " << *itA << std::endl;
          RU._update_barcode(R.get_pivot(*itB), *itA);
          targetInB = false;

          // Target is last element of B, so related bar must be removed
          if ( std::next(itB) == std::prev(B.end()) ) {
            RU.barcode_.erase(indexToBar_[M::get_row_index(*itA)]);
          }

        }

        // Perform addition
        R.add_to(*itB, M::get_row_index(*itA));
        U.add_to(M::get_row_index(*itA), *itB);
        itA--;
      }

      // Clean up
      R.zero_column(columnIndex);
      U.zero_column(columnIndex);
      removedColumns.insert(columnIndex);

      // R.print();
      // U.print();
    }
  };

  // Part of Step 1: Get inverse row
  static Column get_inverse_row(Master_RU_matrix& RU, Index r, std::set<Index>& removedColumns) {
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

    // Get rid of unnecessary entries
    // for (auto i : removedColumns) {
    //   res.clear(i);
    // }

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

    // Swap all in one sweep
    if (removedColumns.size() > 0) {
      auto it_rc = removedColumns.begin();
      bool at_end = false;
      for (int i = *it_rc; i < size; i++) {
        if (!at_end && *it_rc == i) {
          shift += 1;
          if (it_rc == removedColumns.end()) {
            at_end = true;
          } else {
            it_rc++;
          }
        } else {
          U.swap_rows(i - shift, i);
          // std::cout << "*swapping " << i - shift << " and " << i << "\n";
        }
      }
    }

    // Erase all in one sweep
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
  static void _update_persistence(Master_RU_matrix& RU, std::set<Index>& removedColumns){

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
