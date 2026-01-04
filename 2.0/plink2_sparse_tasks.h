// This file is part of PLINK 2.0, copyright (C) 2005-2025 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef __PLINK2_SPARSE_TASKS_H__
#define __PLINK2_SPARSE_TASKS_H__

#include "include/plink2_base.h"

#ifdef __cplusplus
namespace plink2 {
#endif

// Task file entry (20 bytes)
typedef struct TaskEntryStruct {
  uint32_t variant_id;
  uint32_t trait_id;
  double coefficient;
  uint32_t padding;  // For 20-byte alignment
} TaskEntry;

// Task metadata from index file
typedef struct TaskMetadataStruct {
  uint32_t task_id;
  uint64_t num_entries;
  uint32_t trait_start;
  uint32_t trait_end;
  char filename[256];
} TaskMetadata;

// Task index loaded from .index file
typedef struct TaskIndexStruct {
  uint32_t version;
  uint32_t num_tasks;
  uint64_t num_variants;
  uint64_t num_traits;
  TaskMetadata* tasks;  // Array of task metadata
} TaskIndex;

// Task data loaded from .taskN file
typedef struct TaskDataStruct {
  uint32_t task_id;
  uint64_t num_entries;
  uint32_t trait_start;
  uint32_t trait_end;
  TaskEntry* entries;  // Array of task entries
} TaskData;

// Function declarations
PglErr LoadTaskIndex(const char* index_filename, TaskIndex* index);
void FreeTaskIndex(TaskIndex* index);

PglErr LoadTaskData(const char* task_filename, TaskData* data);
void FreeTaskData(TaskData* data);

PglErr ScoreSparseTaskBased(
    const TaskIndex* index,
    const char* tasks_prefix,
    const double* dosages_vmaj,      // [variants × samples]
    double* scores_cmaj,              // [traits × samples]  
    uint32_t num_samples,
    uint32_t num_threads
);

#ifdef __cplusplus
}  // namespace plink2
#endif

#endif  // __PLINK2_SPARSE_TASKS_H__
