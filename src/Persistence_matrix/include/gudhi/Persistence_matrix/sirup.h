#ifndef SIRUP_METHODS_H
#define SIRUP_METHODS_H

#include <gudhi/persistence_matrix_options.h>
#include <gudhi/Persistence_matrix/RU_matrix.h>
#include <gudhi/Persistence_matrix/ru_pairing.h>
#include <gudhi/Persistence_matrix/columns/entry_types.h>

namespace Gudhi {
namespace persistence_matrix {

template <class Master_matrix>
class SiRUP_methods {
  using Column = typename Master_matrix::Column;
  using Index = typename Master_matrix::Index;

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;

 public:
  static void TEST_remove_maximal_cell(Master_RU_matrix& RU, std::vector<Index> removeColumns, bool is_index_updated) {
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;

    std::cout << "*num cols: " << R.matrix_.size() << std::endl;

    auto size = R.get_number_of_columns();
    R.print();
    U.print();

    R.zero_column(removeColumns[0]);
    U.zero_column(removeColumns[0]);
    R.erase_empty_column(removeColumns[0]);
    U.erase_empty_column(removeColumns[0]);

    for (int i = removeColumns[0]; i < size - 1; i++) {
      U.swap_rows(i, i + 1);
      R.swap_rows(i, i + 1);
    }

    std::cout << "--------------------------------\n\n";
    R.print();
    U.print();

    RU.barcode_.clear();
    RU.indexToBar_.clear();

    std::cout << "*num cols: " << R.matrix_.size() << std::endl;

    for (Index i = 0; i < R.get_number_of_columns(); i++) {
      std::cout << "*col " << i << ": ";
      if (!(R.is_zero_column(i))) {
        std::cout << "pivot on row " << R.get_pivot(i) << std::endl;
        RU._update_barcode(R.get_pivot(i), i);
      } else {
        std::cout << "zero col" << std::endl;
        RU._add_bar(RU.get_column_dimension(i), i);
      }
    }
  }

  static void remove_maximal_cell(Master_RU_matrix& RU, std::vector<Index> removeColumns, bool is_index_updated) {
    static_assert(Master_matrix::Option_list::has_removable_columns,
                  "'remove_maximal_cell' is not implemented for the chosen options.");

    // Shorthand and containers
    auto op = RU.operators_;
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;
    using M = Master_matrix;
    constexpr const Index nullDeath = M::template get_null_value<Index>();
    std::set<Index> removedColumns;

    // Step 1: Execute SiRUP
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

      // Get row of V (from its inverse U)
      // Note this differs from "A" in paper, as this A contains diagonal element
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

      // Construct vector of earliest-highest pivot columns
      std::vector<Index> B = {columnIndex};
      auto itA = A.begin();
      itA++;
      while (itA != A.end() && !R.is_zero_column(M::get_row_index(*itA))) {
        Index b = R.get_pivot(M::get_row_index(*itA));
        if (b < R.get_pivot(B.back())) {
          B.push_back(M::get_row_index(*itA));
        }
        itA++;
      }
      if (itA != A.end()) {
        B.push_back(M::get_row_index(*itA));
      }

      // Perform column additions and barcode reindexation
      // If column has not been added to other columns, only update persistence
      if (A.size() == 1) {
        if (!R.is_zero_column(columnIndex)) {
          auto& itBar = RU.indexToBar_.at(columnIndex);
          itBar->death = nullDeath;
          RU.indexToBar_.try_emplace(itBar->birth, itBar);  // Ask Hannah: is this necessary?
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

          // Perform addition
          R.add_to(*itB, M::get_row_index(*itA));
          U.add_to(M::get_row_index(*itA), *itB);
          itA--;
        }
      }

      // Ensure columns are zero
      R.zero_column(columnIndex);
      U.zero_column(columnIndex);
      removedColumns.insert(columnIndex);
    }

    // Step 2: Clear the rows and columns
    auto size = R.get_number_of_columns();
    int shift = 0;

    // Swap all rows in one sweep
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
        } else if ( shift > 0) {
          U.swap_rows(i - shift, i);
          R.swap_rows(i - shift, i);
        }
      }
    }

    // Erase all columns in one sweep
    shift = 0;
    for (Index columnIndex : removedColumns) {
      R.erase_empty_column(columnIndex - shift);
      U.erase_empty_column(columnIndex - shift);
      shift++;
    }

    // Step 3: Recompute barcode
    RU.barcode_.clear();
    RU.indexToBar_.clear();

    for (Index i = 0; i < R.get_number_of_columns() - removedColumns.size(); i++) {
      if (!(R.is_zero_column(i))) {
        auto temp = R.get_pivot(i);

        // ERROR CATCH START
        if (temp > R.get_number_of_columns() - removedColumns.size()) {
          std::cout << "* hi index at i=" << i << " and piv=" << temp << std::endl;
          std::cout << "* expected matrix size is " << R.get_number_of_columns() - removedColumns.size() << " rows and cols" << std::endl;
        }
        // ERROR CATCH END

        auto& barIt = RU.indexToBar_.at(R.get_pivot(i));
        barIt->death = i;
        RU.indexToBar_.try_emplace(i, barIt);
      } else {
        RU._add_bar(RU.get_column_dimension(i), i);
        RU.barcode_.back().death = nullDeath;
      }
    }

    // DEBUG START
    // std::cout << "------------------------------------------\n";
    // R.print();
    // U.print();
    // DEBUG END

  };
};

}  // namespace persistence_matrix
}  // namespace Gudhi

#endif  // SIRUP_METHODS_H
