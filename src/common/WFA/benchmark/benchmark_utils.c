/*
 *                             The MIT License
 *
 * Wavefront Alignments Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignments Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignments Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: Benchmark utils
 */

#include "benchmark/benchmark_utils.h"
#include "edit/edit_table.h"
#include "edit/edit_dp.h"
#include "gap_lineal/nw.h"
#include "gap_affine/affine_table.h"
#include "gap_affine/swg.h"

/*
 * Setup
 */
void benchmark_align_input_clear(
    align_input_t* const align_input) {
  // Accuracy Stats
  counter_reset(&(align_input->align));
  counter_reset(&(align_input->align_correct));
  counter_reset(&(align_input->align_score));
  counter_reset(&(align_input->align_score_total));
  counter_reset(&(align_input->align_score_diff));
  counter_reset(&(align_input->align_cigar));
  counter_reset(&(align_input->align_bases));
  counter_reset(&(align_input->align_matches));
  counter_reset(&(align_input->align_mismatches));
  counter_reset(&(align_input->align_del));
  counter_reset(&(align_input->align_ins));
  wavefronts_stats_clear(&(align_input->wavefronts_stats));
}
/*
 * Check
 */
void benchmark_check_alignment(
    align_input_t* const align_input,
    edit_cigar_t* const edit_cigar_computed) {
  // Compute correct CIGAR
  if ((align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) ||
      (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT)) {
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_AFFINE) {
      // Compute correct
      affine_table_t affine_table;
      affine_table_allocate(
          &affine_table,align_input->pattern_length,
          align_input->text_length,align_input->mm_allocator);
      if (align_input->check_bandwidth < 0) {
        swg_compute(&affine_table,align_input->check_affine_penalties,
            align_input->pattern,align_input->pattern_length,
            align_input->text,align_input->text_length);
      } else {
        swg_compute_banded(&affine_table,align_input->check_affine_penalties,
            align_input->pattern,align_input->pattern_length,
            align_input->text,align_input->text_length,
            align_input->check_bandwidth);
      }
      const int score_correct = edit_cigar_score_gap_affine(
          &affine_table.edit_cigar,align_input->check_affine_penalties);
      const int score_computed = edit_cigar_score_gap_affine(
          edit_cigar_computed,align_input->check_affine_penalties);
      // Check alignment
      benchmark_check_alignment_using_template(
          align_input,edit_cigar_computed,score_computed,
          &affine_table.edit_cigar,score_correct);
      // Free
      affine_table_free(&affine_table,align_input->mm_allocator);
    } else if(align_input->debug_flags & ALIGN_DEBUG_CHECK_DISTANCE_METRIC_GAP_LINEAL) {
      // Compute correct
      edit_table_t edit_table;
      edit_table_allocate(
          &edit_table,align_input->pattern_length,
          align_input->text_length,align_input->mm_allocator);
      nw_compute(&edit_table,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length,align_input->check_lineal_penalties);
      const int score_correct = edit_cigar_score_gap_lineal(
          &edit_table.edit_cigar,align_input->check_lineal_penalties);
      const int score_computed = edit_cigar_score_gap_lineal(
          edit_cigar_computed,align_input->check_lineal_penalties);
      // Check alignment
      benchmark_check_alignment_using_template(
          align_input,edit_cigar_computed,score_computed,
          &edit_table.edit_cigar,score_correct);
      // Free
      edit_table_free(&edit_table,align_input->mm_allocator);
    } else { // ALIGN_DEBUG_CHECK_DISTANCE_METRIC_EDIT
      edit_table_t edit_table;
      edit_table_allocate(
          &edit_table,align_input->pattern_length,
          align_input->text_length,align_input->mm_allocator);
      if (align_input->check_bandwidth < 0) {
        edit_dp_compute(&edit_table,
            align_input->pattern,align_input->pattern_length,
            align_input->text,align_input->text_length);
      } else {
        edit_dp_compute_banded(&edit_table,
            align_input->pattern,align_input->pattern_length,
            align_input->text,align_input->text_length,
            align_input->check_bandwidth);
      }
      const int score_correct = edit_cigar_score_edit(&edit_table.edit_cigar);
      const int score_computed = edit_cigar_score_edit(edit_cigar_computed);
      // Check alignment
      benchmark_check_alignment_using_template(
          align_input,edit_cigar_computed,score_computed,
          &edit_table.edit_cigar,score_correct);
      // Free
      edit_table_free(&edit_table,align_input->mm_allocator);
    }
  } else {
    // Delegate check alignment
    const int score_computed = edit_cigar_score_edit(edit_cigar_computed);
    benchmark_check_alignment_using_template(
        align_input,edit_cigar_computed,score_computed,NULL,-1);
  }
}
void benchmark_check_alignment_using_template(
    align_input_t* const align_input,
    edit_cigar_t* const edit_cigar_computed,
    const int score_computed,
    edit_cigar_t* const edit_cigar_correct,
    const int score_correct) {
  counter_add(&(align_input->align),1);
  counter_add(&(align_input->align_score_total),ABS(score_computed));
  // Debug
  if (align_input->debug_flags) {
    // Display info
    if (align_input->debug_flags & ALIGN_DEBUG_DISPLAY_INFO) {
      benchmark_print_alignment(stderr,align_input,score_computed,edit_cigar_computed,-1,NULL);
    }
    // Check correct
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_CORRECT) {
      bool correct = edit_cigar_check_alignment(stderr,
          align_input->pattern,align_input->pattern_length,
          align_input->text,align_input->text_length,
          edit_cigar_computed,align_input->verbose);
      if (!correct) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INCORRECT ALIGNMENT\n");
          benchmark_print_alignment(stderr,align_input,-1,edit_cigar_computed,-1,NULL);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_correct),1);
      }
      // CIGAR Stats
      int i;
      counter_add(&(align_input->align_bases),align_input->pattern_length);
      for (i=edit_cigar_computed->begin_offset;i<edit_cigar_computed->end_offset;++i) {
        switch (edit_cigar_computed->operations[i]) {
          case 'M': counter_add(&(align_input->align_matches),1); break;
          case 'X': counter_add(&(align_input->align_mismatches),1); break;
          case 'I': counter_add(&(align_input->align_ins),1); break;
          case 'D': default: counter_add(&(align_input->align_del),1); break;
        }
      }
    }
    // Check score
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_SCORE) {
      if (score_computed != score_correct) {
        // Print
        if (align_input->verbose) {
          benchmark_print_alignment(
              stderr,align_input,
              score_computed,edit_cigar_computed,
              score_correct,edit_cigar_correct);
          fprintf(stderr,"(#%d)\t INACCURATE SCORE computed=%d\tcorrect=%d\n",
              align_input->sequence_id,score_computed,score_correct);
        }
        counter_add(&(align_input->align_score_diff),ABS(score_computed-score_correct));
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_score),1);
      }
    }
    // Check alignment
    if (align_input->debug_flags & ALIGN_DEBUG_CHECK_ALIGNMENT) {
      if (edit_cigar_cmp(edit_cigar_computed,edit_cigar_correct) != 0) {
        // Print
        if (align_input->verbose) {
          fprintf(stderr,"INACCURATE ALIGNMENT\n");
          benchmark_print_alignment(
              stderr,align_input,
              -1,edit_cigar_computed,
              -1,edit_cigar_correct);
        }
        // Quit
        return;
      } else {
        counter_add(&(align_input->align_cigar),1);
      }
    }
  }
}
/*
 * Display
 */
void benchmark_print_alignment(
    FILE* const stream,
    align_input_t* const align_input,
    const int score_computed,
    edit_cigar_t* const cigar_computed,
    const int score_correct,
    edit_cigar_t* const cigar_correct) {
  // Print Sequence
  fprintf(stream,"ALIGNMENT (#%d)\n",align_input->sequence_id);
  fprintf(stream,"  PATTERN  %s\n",align_input->pattern);
  fprintf(stream,"  TEXT     %s\n",align_input->text);
  // Print CIGARS
  if (cigar_computed != NULL && score_computed != -1) {
    fprintf(stream,"    COMPUTED\tscore=%d\t",score_computed);
    edit_cigar_print(stream,cigar_computed);
    fprintf(stream,"\n");
  }
  if (cigar_computed != NULL) {
    edit_cigar_print_pretty(stream,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        cigar_computed,align_input->mm_allocator);
  }
  if (cigar_correct != NULL && score_correct != -1) {
    fprintf(stream,"    CORRECT \tscore=%d\t",score_correct);
    edit_cigar_print(stream,cigar_correct);
    fprintf(stream,"\n");
  }
  if (cigar_correct != NULL) {
    edit_cigar_print_pretty(stream,
        align_input->pattern,align_input->pattern_length,
        align_input->text,align_input->text_length,
        cigar_correct,align_input->mm_allocator);
  }
}
/*
 * Stats
 */
void benchmark_print_stats(
    FILE* const stream,
    align_input_t* const align_input,
    const bool print_wf_stats) {
  // General stats
  fprintf(stream,"[Accuracy]\n");
  fprintf(stream," => Alignments.Correct     ");
  counter_print(stream,&align_input->align_correct,&align_input->align,"alg       ",true);
  fprintf(stream," => Score.Correct          ");
  counter_print(stream,&align_input->align_score,&align_input->align,"alg       ",true);
  fprintf(stream,"   => Score.Total          ");
  counter_print(stream,&align_input->align_score_total,NULL,"score uds.",true);
  fprintf(stream,"     => Score.Diff         ");
  counter_print(stream,&align_input->align_score_diff,&align_input->align_score_total,"score uds.",true);
  fprintf(stream," => CIGAR.Correct          ");
  counter_print(stream,&align_input->align_cigar,&align_input->align,"alg       ",true);
  fprintf(stream,"   => CIGAR.Matches        ");
  counter_print(stream,&align_input->align_matches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Mismatches     ");
  counter_print(stream,&align_input->align_mismatches,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Insertions     ");
  counter_print(stream,&align_input->align_ins,&align_input->align_bases,"bases     ",true);
  fprintf(stream,"   => CIGAR.Deletions      ");
  counter_print(stream,&align_input->align_del,&align_input->align_bases,"bases     ",true);
}






