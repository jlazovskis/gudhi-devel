#ifndef SIRUP_METHODS_H
#define SIRUP_METHODS_H

#include <gudhi/Persistence_matrix/RU_matrix.h>
#include <gudhi/Persistence_matrix/columns/entry_types.h>
#include <gudhi/Persistence_matrix/ru_pairing.h>
#include <gudhi/persistence_matrix_options.h>

namespace Gudhi {
namespace persistence_matrix {

template <class Master_matrix>
class SiRUP_methods {
  using Column = typename Master_matrix::Column;
  using Index = typename Master_matrix::Index;
  using Field_element = typename Master_matrix::Element;

 private:
  using Master_RU_matrix = typename Master_matrix::Master_RU_matrix;
  class binary_indexed_tree {
   private:
    Index n;
    std::vector<Index> tree;

   public:
    binary_indexed_tree(Index n) : n(n), tree(n + 1, 0) {
      for (Index i = 1; i <= n; ++i) {
        add(i, 1);
      }
    }
    void add(Index i, Index diff) {
      while (i <= n) {
        tree[i] += diff;
        i += i & -i;
      }
    }
    void remove(Index k) { add(k + 1, -1); }
    Index find(Index k) {
      Index pos = 0;
      for (Index bit = 1 << 30; bit > 0; bit >>= 1) {
        Index next = pos + bit;
        if (next <= n && tree[next] <= k) {
          pos = next;
          k -= tree[next];
        }
      }
      return pos;
    }
  };
  // static Field_element element_at_index(Column col, Index idx) {
  //   using M = Master_matrix;
  //   bool found = false;
  //   auto it = col.begin();
  //   while (!found && it!=col.end()){
  //       if ( M::get_row_index(*it) == idx ) {
  //           found = true;
  //       }
  //       it++;
  //   }
  //   if (found) {
  //       return M::get_element(*std::prev(it));
  //   } else {
  //       return (Field_element) 0;
  //   }
  // };

 public:
  static void remove_maximal_cell(Master_RU_matrix& RU, std::vector<Index> removeColumns, bool is_index_updated) {
    static_assert(Master_matrix::Option_list::has_removable_columns,
                  "'remove_maximal_cell' is not implemented for the chosen options.");

    // Shorthand and containers
    auto op = RU.operators_;
    auto R = RU.reducedMatrixR_;
    auto U = RU.mirrorMatrixU_;
    using M = Master_matrix;
    std::set<Index> removedColumns;

    // Constant values
    constexpr const Index nullDeath = M::template get_null_value<Index>();
    const auto size_input = R.get_number_of_columns();
    const auto size_output = R.get_number_of_columns() - removeColumns.size();

    // DEBUG START
    // std::cout << "* removeColumns" << std::endl;
    // for (auto it_rc = removeColumns.begin(); it_rc != removeColumns.end();
    //     ++it_rc) {
    //     std::cout << "* i=" << *it_rc << std::endl;
    // }
    // DEBUG END

    // Step 0: Reindex and reorder
    // Ensure indexing is absolute
    std::vector<Index> removeColumns_ordered;
    if (is_index_updated) {
      Index maxindex = removeColumns[0];
      for (int i = 1; i < removeColumns.size(); i++) {
        if (removeColumns[i] > maxindex) {
          maxindex = removeColumns[i];
        }
      }
      binary_indexed_tree bit = binary_indexed_tree(maxindex + removeColumns.size() + 1);
      for (Index k : removeColumns) {
        Index k_original = bit.find(k);
        removeColumns_ordered.push_back(k_original);
        bit.remove(k_original);
      }
    } else {
      removeColumns_ordered = removeColumns;
    };

    // Ensure highest dimensions first, earliest indices first
    std::sort(removeColumns_ordered.begin(), removeColumns_ordered.end(), [&RU](const Index& a, const Index& b) {
      if (RU.get_column_dimension(a) == RU.get_column_dimension(b)) {
        return a < b;
      } else {
        return RU.get_column_dimension(a) > RU.get_column_dimension(b);
      }
    });

    // DEBUG START
    // std::cout << "* removedColumns_ordered" << std::endl;
    // for (auto it_rc = removeColumns_ordered.begin(); it_rc != removeColumns_ordered.end(); ++it_rc) {
    //   std::cout << *it_rc << " (" << RU.get_column_dimension(*it_rc) << ") ";
    //   // std::cout << *it_rc << std::endl;
    // }
    // std::cout << std::endl;
    // DEBUG END

    // DEBUG START
    // for (Index i = 0; i < size_input; i++) {
    //   if (!(R.is_zero_column(i))) {
    //     auto temp = R.get_pivot(i);
    //     if (temp > size_output) {
    //       std::cout << "** BOOPS: before sirup col " << i << " has pivot " << temp << " (>" << size_output << ")"
    //                 << std::endl;
    //     }
    //   }
    // }
    // DEBUG END

    // DEBUG START
    // for (Index i = 0; i < size_input; i++) {
    //   const Column col = Column(R.get_column(i), nullptr);
    //   for ( auto it=col.begin(); it!=col.end(); it++) {
    //     bool found = false;
    //      for ( auto it2=removeColumns_ordered.begin(); it2 !=removeColumns_ordered.end(); it2++) {
    //         if ( M::get_row_index(*it) == *it2) {
    //             found = true;
    //         }
    //     }
    //     if ( found ) { std::cout << "** Column " << i << " (before) has nonzero entry in row-to-be-removed" <<
    //     std::endl;}
    //   }
    // }
    // DEBUG END

    // Step 1: Execute SiRUP
    // Iterate over columns requested to remove
    for (int i = 0; i < removeColumns_ordered.size(); i++) {
      Index columnIndex = removeColumns_ordered[i];
      auto columnDim = RU.get_column_dimension(columnIndex);

      // DEBUG START
      // std::cout << " * row 3 (before adding " << columnIndex << ") is nonzero in col: ";
      // for (int j = 0; j < size_input; j++) {
      //   const Column col = Column(R.get_column(j), nullptr);
      //   for (auto it = col.begin(); it != col.end(); it++) {
      //     if (M::get_row_index(*it) == 3) {
      //       std::cout << j << " ";
      //     }
      //   }
      // }
      // std::cout << std::endl;
      // DEBUG END

      // DEBUG START
      // std::cout << " * columnIndex=" << columnIndex << std::endl;
      // DEBUG END

      // DEBUG START
      // std::cout << " ** row of U at index " << columnIndex << " is ";
      // for ( auto itt = A_const.begin(); itt!=A_const.end(); itt++){
      //   std::cout << M::get_row_index(*itt) << " ";
      // }
      // std::cout << std::endl;
      // DEBUG END

      // Get row of V (from its inverse U)
      // Note this differs from "A" in paper, as this A contains diagonal
      // element
      Column A = Column(U.get_column(columnIndex), nullptr);
      std::vector<Index> A_index = {columnIndex};
      std::vector<Field_element> A_entry = {1};
      for (Index col = columnIndex + 1; col < size_input; ++col) {
          auto itA = A.begin();
          while (itA != A.end() && M::get_row_index(*itA) < col) { ++itA; }
          if (itA != A.end() && M::get_row_index(*itA) == col) {
            // Record current values
            Field_element entry = op->get_characteristic() - M::get_element(*itA);
            A_index.push_back(col);
            A_entry.push_back(entry);
            // Add current row
            Column cur_row = Column(U.get_column(col), nullptr);
            cur_row *= entry;
            A += cur_row;
          }
      }

      // DEBUG START
      // std::cout << " ** at U_row_index=" << U_row_index << ", row " << columnIndex << " of V is ";
      // auto itt = A.begin();
      // for (int intt=columnIndex; intt<size_input; intt++){
      //     if ( itt!=A.end() && intt == M::get_row_index(*itt) ){
      //         std::cout << "1 ";
      //         itt++;
      //     } else {
      //         std::cout << "0 ";
      //     }
      // }
      // std::cout << std::endl;
      // DEBUG END

      // auto itAc = std::next(A_const.begin());
      // while (itAc != A_const.end()) {
      //   Column cur_row = Column(U.get_column(*itAc), nullptr);
      //   cur_row.clear(*itAc);
      //   cur_row *= op->get_characteristic() - M::get_element(*itAc);
      //   A += cur_row;
      //   ++itAc;
      // }

      // DEBUG START
      // if (columnIndex == 4) {
      //   std::cout << " * A for i=4: ";
      //   for (auto it = A.begin(); it != A.end(); ++it) {
      //     std::cout << *it << " (" << RU.get_column_dimension(*it) << ") ";
      //   }
      //   std::cout << std::endl;
      // }
      // DEBUG END

      // Construct vector of earliest-highest pivot columns
      std::vector<Index> B = {columnIndex};
      auto iitA = A_index.begin();
      iitA++;
      while ( iitA!=A_index.end() && !R.is_zero_column(*iitA) ) {
        Index b = R.get_pivot(*iitA);
        if (b < R.get_pivot(B.back())) {
          B.push_back(*iitA);
        }
        iitA++;
      }
      if (iitA != A_index.end()) {
        B.push_back(*iitA);
      }

      // DEBUG START
      // if (columnIndex == 570) {
      //   std::cout << " * B for i=570: ";
      //   for (auto it = B.begin(); it != B.end(); ++it) {
      //     std::cout << *it << " (" << RU.get_column_dimension(*it) << ") ";
      //   }
      //   std::cout << std::endl;
      // }
      // DEBUG END

      // Perform column additions and barcode reindexation
      // Case 1: column has not been added to other columns. Only update persistence.
      if (A_index.size() == 1) {
        if (!R.is_zero_column(columnIndex)) {
          auto& itBar = RU.indexToBar_.at(columnIndex);
          itBar->death = nullDeath;
          RU.indexToBar_.try_emplace(itBar->birth,
                                     itBar);  // Ask Hannah: is this necessary?
        } else {
          RU.barcode_.erase(RU.indexToBar_[columnIndex]);
        }
      }
      // Case 2: column has been added to other columns. Perform column addition.
      else {
        iitA = A_index.end();
        iitA--;
        bool targetInB = false;
        bool targetIsZeroCol = false;
        while (iitA != A_index.begin()) {
          auto itB = B.begin();

          // Case 2.1: target is zero column
          if (R.is_zero_column(*iitA)) {
            itB = std::prev(B.end());
            if (*itB == *iitA) {
              itB--;
              targetInB = true;
            }
            targetIsZeroCol = true;
          }

          // Case 2.2: target is nonzero column
          else {
            Index k = R.get_pivot(*iitA);
            while ((R.get_pivot(*iitA) > k) && (R.get_pivot(*itB) != -1)) {
              itB++;
            }
            if (*itB == M::get_row_index(*iitA)) {
              itB--;
              targetInB = true;
            }
          }

          // DEBUG START
          // std::cout << " * addition for col " << columnIndex << ": adding " << *itB << " to " <<
          // M::get_row_index(*iitA) << std::endl; DEBUG END

          // Perform addition (backwards on U because same matrix on left)
          R.add_to(*itB, *iitA);
          U.add_to(*iitA, *itB);
          iitA--;
        }
      }

      // Ensure columns are zero
      R.zero_column(columnIndex);
      U.zero_column(columnIndex);
      removedColumns.insert(columnIndex);

      // DEBUG START
      //   std::cout << " * col " << columnIndex << ": ";
      //   for ( int j=0; j<size_input; j++) {
      //     const Column col = Column(R.get_column(j), nullptr);
      //     for ( auto it=col.begin(); it!=col.end(); it++) {
      //         if ( M::get_row_index(*it) == columnIndex ){
      //           std::cout << j << " ";
      //         }
      //     }
      // }
      //   std::cout << std::endl;
      // DEBUG END
    }

    // DEBUG START
    // std::cout << "* removedColumns" << std::endl;
    // for (auto it_rc = removedColumns.begin(); it_rc != removedColumns.end();
    //     ++it_rc) {
    //     std::cout << *it_rc << " (" << RU.get_column_dimension(*it_rc) << ") ";
    // }
    // std::cout << std::endl;
    // DEBUG END

    // DEBUG START
    // for (Index i = 0; i < size_input; i++) {
    //   const Column col = Column(R.get_column(i), nullptr);
    //   for (auto it = col.begin(); it != col.end(); it++) {
    //     bool found = false;
    //     Index foundIndex;
    //     for (auto it_rc = removedColumns.begin(); it_rc != removedColumns.end(); it_rc++) {
    //       if (*it_rc == M::get_row_index(*it)) {
    //         found = true;
    //         foundIndex = M::get_row_index(*it);
    //       }
    //     }
    //     if (found) {
    //       std::cout << "** BOOPS: before row swap col " << i << " has nonzero entry in row " << foundIndex <<
    //       std::endl;
    //     }
    //   }
    // }
    // std::cout << " --------------------------------\nSWAPPING" << std::endl;
    // DEBUG END

    // Step 2: Clear the rows and columns
    int shift = 0;

    // Swap all rows in one sweep
    if (removedColumns.size() > 0) {
      auto it_rc = removedColumns.begin();
      bool at_rc_end = false;
      for (int i = *it_rc; i < size_input; i++) {
        if (!at_rc_end && *it_rc == i) {
          shift += 1;
          it_rc++;
          if (it_rc == removedColumns.end()) {
            at_rc_end = true;
          }
        } else {
          U.swap_rows(i - shift, i);
          R.swap_rows(i - shift, i);

          // DEBUG START
          // std::cout << i-shift << " <-> " << i << " ";
          // DEBUG END
        }
      }
    }
    // std::cout << std::endl;

    // DEBUG START
    // for (Index i = 0; i < size_input; i++) {
    //   if (!(R.is_zero_column(i))) {
    //     auto temp = R.get_pivot(i);
    //     if (temp > size_output) {
    //       std::cout << "** OOPS: after row swap col " << i << " (" << RU.get_column_dimension(i) << ") has pivot "
    //                 << temp << " (" << RU.get_column_dimension(temp) << ") (>" << size_output << ")" << std::endl;
    //     }
    //   }
    // }
    // DEBUG END

    // Erase all columns in one sweep
    shift = 0;
    for (Index columnIndex : removedColumns) {
      R.erase_empty_column(columnIndex - shift);
      U.erase_empty_column(columnIndex - shift);
      shift++;
    }

    // DEBUG START
    std::cout << "* From " << size_input << "x" << size_input << " to " << size_output << "x" << size_output
              << std::endl;
    // DEBUG END

    // Step 3: Recompute barcode
    RU.barcode_.clear();
    RU.indexToBar_.clear();

    for (Index i = 0; i < size_output; i++) {
      if (!(R.is_zero_column(i))) {
        auto temp = R.get_pivot(i);

        // DEBUG START
        // std::cout << "* i=" << i << " and piv=" << temp << std::endl;
        // DEBUG END

        // ERROR CATCH START
        if (temp > size_output) {
          std::cout << "* hi index at i=" << i << " and piv=" << temp << std::endl;
          std::cout << "* expected matrix size is " << size_output << " rows and cols" << std::endl;
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
