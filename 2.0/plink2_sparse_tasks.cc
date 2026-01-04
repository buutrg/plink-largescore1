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

#include "plink2_sparse_tasks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef __cplusplus
namespace plink2 {
#endif

PglErr LoadTaskIndex(const char* index_filename, TaskIndex* index) {
  FILE* f = fopen(index_filename, "rb");
  if (!f) {
    return kPglRetOpenFail;
  }
  
  // Read header (32 bytes)
  char magic[8];
  if (fread(magic, 1, 8, f) != 8 || memcmp(magic, "PLNKTASK", 8) != 0) {
    fclose(f);
    return kPglRetMalformedInput;
  }
  
  fread(&index->version, sizeof(uint32_t), 1, f);
  fread(&index->num_tasks, sizeof(uint32_t), 1, f);
  fread(&index->num_variants, sizeof(uint64_t), 1, f);
  fread(&index->num_traits, sizeof(uint64_t), 1, f);
  
  // Allocate task metadata array
  index->tasks = (TaskMetadata*)malloc(index->num_tasks * sizeof(TaskMetadata));
  if (!index->tasks) {
    fclose(f);
    return kPglRetNomem;
  }
  
  // Read task metadata (32 bytes each)
  for (uint32_t i = 0; i < index->num_tasks; ++i) {
    TaskMetadata* meta = &index->tasks[i];
    fread(&meta->task_id, sizeof(uint32_t), 1, f);
    fread(&meta->num_entries, sizeof(uint64_t), 1, f);
    fread(&meta->trait_start, sizeof(uint32_t), 1, f);
    fread(&meta->trait_end, sizeof(uint32_t), 1, f);
    
    uint64_t reserved;
    fread(&reserved, sizeof(uint64_t), 1, f);
    
    // Construct filename from task_id (will be set by caller)
    meta->filename[0] = '\0';
  }
  
  fclose(f);
  return kPglRetSuccess;
}

void FreeTaskIndex(TaskIndex* index) {
  if (index->tasks) {
    free(index->tasks);
    index->tasks = nullptr;
  }
}

PglErr LoadTaskData(const char* task_filename, TaskData* data) {
  FILE* f = fopen(task_filename, "rb");
  if (!f) {
    return kPglRetOpenFail;
  }
  
  // Read task header (24 bytes)
  fread(&data->num_entries, sizeof(uint64_t), 1, f);
  fread(&data->trait_start, sizeof(uint32_t), 1, f);
  fread(&data->trait_end, sizeof(uint32_t), 1, f);
  
  uint64_t reserved;
  fread(&reserved, sizeof(uint64_t), 1, f);
  
  // Allocate entries array
  data->entries = (TaskEntry*)malloc(data->num_entries * sizeof(TaskEntry));
  if (!data->entries) {
    fclose(f);
    return kPglRetNomem;
  }
  
  // Read all entries (20 bytes each)
  size_t bytes_to_read = data->num_entries * sizeof(TaskEntry);
  size_t bytes_read = fread(data->entries, 1, bytes_to_read, f);
  
  fclose(f);
  
  if (bytes_read != bytes_to_read) {
    free(data->entries);
    data->entries = nullptr;
    return kPglRetReadFail;
  }
  
  return kPglRetSuccess;
}

void FreeTaskData(TaskData* data) {
  if (data->entries) {
    free(data->entries);
    data->entries = nullptr;
  }
}

// Process a single task: load from disk and accumulate scores
static void ProcessTask(
    const char* task_filename,
    const double* dosages_vmaj,
    double* scores_cmaj,
    uint32_t num_samples
) {
  TaskData data;
  PglErr reterr = LoadTaskData(task_filename, &data);
  if (reterr != kPglRetSuccess) {
    return;  // TODO: Better error handling
  }
  
  // Process all entries for this task
  for (uint64_t i = 0; i < data.num_entries; ++i) {
    uint32_t variant_id = data.entries[i].variant_id;
    uint32_t trait_id = data.entries[i].trait_id;
    double coef = data.entries[i].coefficient;
    
    const double* variant_dosages = &dosages_vmaj[variant_id * num_samples];
    double* trait_scores = &scores_cmaj[trait_id * num_samples];
    
    // Vectorizable accumulation loop
    for (uint32_t s = 0; s < num_samples; ++s) {
      trait_scores[s] += coef * variant_dosages[s];
    }
  }
  
  FreeTaskData(&data);
}

PglErr ScoreSparseTaskBased(
    const TaskIndex* index,
    const char* tasks_prefix,
    const double* dosages_vmaj,
    double* scores_cmaj,
    uint32_t num_samples,
    uint32_t num_threads
) {
  // Initialize scores to zero
  memset(scores_cmaj, 0, index->num_traits * num_samples * sizeof(double));
  
  // OpenMP task-based parallelism
  #pragma omp parallel num_threads(num_threads)
  {
    #pragma omp single
    {
      // Spawn tasks dynamically
      for (uint32_t task_id = 0; task_id < index->num_tasks; ++task_id) {
        #pragma omp task firstprivate(task_id)
        {
          char task_filename[512];
          snprintf(task_filename, sizeof(task_filename), 
                   "%s.task%u", tasks_prefix, task_id);
          
          ProcessTask(task_filename, dosages_vmaj, scores_cmaj, num_samples);
        }
      }
      // Implicit taskwait here (all tasks complete before exiting single)
    }
  }
  
  return kPglRetSuccess;
}

#ifdef __cplusplus
}  // namespace plink2
#endif
