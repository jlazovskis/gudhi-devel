#ifndef SIRUP_METHODS_H
#define SIRUP_METHODS_H

#include <gudhi/Persistence_matrix/RU_matrix.h>
#include <gudhi/persistence_matrix_options.h>
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

 public:
  static void remove_maximal_cell(const typename Master_matrix::Field_operators* op,
                                  typename Master_matrix::Master_boundary_matrix& R,
                                  typename Master_matrix::Master_base_matrix& U, std::vector<Index> removeColumns,
                                  bool is_index_updated) {
    static_assert(Master_matrix::Option_list::has_removable_columns,
                  "'remove_maximal_cell' is not implemented for the chosen options.");

    // Container for columns already removed. Indexation does not change.
    std::set<Index> removedColumns;

    // Steps of algorithm implementation
    _execute_sirup(op, R, U, removeColumns, removedColumns, is_index_updated);
    _clear_rows_columns(op, R, U, removedColumns, is_index_updated);
    _update_barcode(op, R, U, removedColumns);
  }

  // Step 1: Execute SiRUP
  static void _execute_sirup(const typename Master_matrix::Field_operators* op,
                             typename Master_matrix::Master_boundary_matrix& R,
                             typename Master_matrix::Master_base_matrix& U, std::vector<Index> removeColumns,
                             std::set<Index> removedColumns, bool is_index_updated) {
    // Shorthand
    using M = Master_matrix;

    // Iterate over columns requested to remove
    for (int i = 0; i < removeColumns.size(); i++) {
      Index columnIndex = removeColumns[i];
      if (is_index_updated && removedColumns.size() > 0) {
        for (auto it_rc = removedColumns.begin(); it_rc != removedColumns.end(); ++it_rc){
          if ( *it_rc <= columnIndex ) {
            columnIndex += 1;
          }
        }
      }
      std::cout << "removing at index " << columnIndex << std::endl;

      // Definitions
      Index size = R.get_number_of_columns();
      // std::cout << "* size: " << size <<"\n";

      // Get row of V (by inverting U)
      const auto A = get_inverse_row(op, U, columnIndex, removedColumns);
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
      while (itA != A.begin()) {
        auto itB = B.begin();
        // std::cout << "*considering col " << M::get_row_index(*itA);
        // if (R.is_zero_column(M::get_row_index(*itA))) { std:: cout << " (is zero)\n";}
        // else { std:: cout << " (is not zero)\n"; }

        // Separate case for zero columns
        if (R.is_zero_column(M::get_row_index(*itA))) {
          // if (R.is_zero_column(row_shift(M::get_row_index(*itA)))) {
          itB = std::prev(B.end());
          if (*itB == M::get_row_index(*itA)) {
            // if (*itB == row_shift(M::get_row_index(*itA))) {
            itB--;
          }
        }

        // Separate case for nonzero columns
        else {
          Index k = R.get_pivot(M::get_row_index(*itA));
          // Index k = R.get_pivot(row_shift(M::get_row_index(*itA)));
          while ((R.get_pivot(*itB) > k) && (R.get_pivot(*itB) != -1)) {
            itB++;
          }
          if (*itB == M::get_row_index(*itA)) {
            // if (*itB == row_shift(M::get_row_index(*itA))) {
            itB--;
          }
        }

        // Perform addition and go to next step
        // std::cout << "adding " << *itB << " to " << M::get_row_index(*itA) << std::endl;
        // std::cout << "adding " << *itB << " to " << row_shift(M::get_row_index(*itA)) << std::endl;
        R.add_to(*itB, M::get_row_index(*itA));
        // R.add_to(*itB, row_shift(M::get_row_index(*itA)));
        U.add_to(M::get_row_index(*itA), *itB);
        // U.add_to(*itB, row_shift(M::get_row_index(*itA)));
        itA--;
      }

      // Clean up
      R.zero_column(columnIndex);
      U.zero_column(columnIndex);
      removedColumns.insert(columnIndex);
    }
  };

  // Part of Step 1: Get inverse row
  static Column get_inverse_row(const typename Master_matrix::Field_operators* op,
                                typename Master_matrix::Master_base_matrix& U, Index r,
                                std::set<Index> removedColumns) {
    // Shorthand
    using M = Master_matrix;
    using E = typename M::Element;

    // Definitions
    Column res = Column(U.get_column(r), nullptr);
    auto size = U.get_number_of_columns();
    auto it = std::next(res.begin());
    auto it_rc = removedColumns.begin();
    if (removedColumns.size() > 0) {
      it_rc++;
    }

    // Construct inverse
    while (it != res.end()) {
      if ((it_rc != removedColumns.end()) && (*it_rc == M::get_element(*it))) {
        ++it_rc;
      } else {
        Column cur_row = Column(U.get_column(*it), nullptr);
        cur_row.clear(*it);
        cur_row *= op->get_characteristic() - M::get_element(*it);
        res += cur_row;
      }
      ++it;
    }

    // Get rid of unnecessary entries
    for (auto i : removedColumns) {
      res.clear(i);
    }

    // Return result
    return res;
  };

  // Step 2: Clear the rows and columns
  static void _clear_rows_columns(const typename Master_matrix::Field_operators* op,
                                  typename Master_matrix::Master_boundary_matrix& R,
                                  typename Master_matrix::Master_base_matrix& U, std::set<Index> removedColumns,
                                  bool is_index_updated) {
    // Remove columns
    // R.erase_empty_column(columnIndex);
    // U.erase_empty_column(columnIndex);

    auto size = R.get_number_of_columns();
    int shift = 0;

    // Swap all in one sweep
    if (removedColumns.size() > 0) {
      auto it_rc = std::next(removedColumns.begin());
      bool at_end = false;
      for (int i = *it_rc; i < size - 1; i++) {
        if (!at_end && *it_rc == i) {
          shift += 1;
          if (it_rc == removedColumns.end()) {
            at_end = true;
          } else {
            it_rc++;
          }
        } else {
          U.swap_rows(i - shift, i);
        }

        // Erase all in one sweep
        for (Index columnIndex : removedColumns) {
          R.erase_empty_column(columnIndex);
          U.erase_empty_column(columnIndex);
          U.erase_empty_row(size - 1);
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
      }
    };
  };

  // Step 3: Update the barcode
  static void _update_barcode(const typename Master_matrix::Field_operators* op,
                              typename Master_matrix::Master_boundary_matrix& R,
                              typename Master_matrix::Master_base_matrix& U, std::set<Index> removedColumns){};
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // SIRUP_METHODS_H
