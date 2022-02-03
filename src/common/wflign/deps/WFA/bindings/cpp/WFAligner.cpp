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
 * DESCRIPTION: C++ bindings for the WaveFront Alignment modules
 */

#include "WFAligner.hpp"

extern "C" {
  #include "../../wavefront/wavefront_align.h"
}

/*
 * General Wavefront Aligner
 */
WFAligner::WFAligner(
    const bool onlyScore,
    const MemoryModel memoryModel) {
  this->attributes = wavefront_aligner_attr_default;
  switch (memoryModel) {
    case WavefrontMemoryHigh: this->attributes.memory_mode = wavefront_memory_high; break;
    case WavefrontMemoryMed: this->attributes.memory_mode = wavefront_memory_med; break;
    case WavefrontMemoryLow: this->attributes.memory_mode = wavefront_memory_low; break;
    default: this->attributes.memory_mode = wavefront_memory_full; break;
  }
  this->attributes.alignment_scope = onlyScore ? compute_score : compute_alignment;
  this->wfAligner = nullptr;
}
WFAligner::~WFAligner() {
  wavefront_aligner_delete(wfAligner);
}
/*
 * Configuration
 */
void WFAligner::setReductionNone() {
  wavefront_aligner_set_reduction_none(wfAligner);
}
void WFAligner::setReductionAdaptive(
    const int min_wavefront_length,
    const int max_distance_threshold) {
  wavefront_aligner_set_reduction_adaptive(
      wfAligner,min_wavefront_length,max_distance_threshold);
}
void WFAligner::setMatchFunct(
    int (*matchFunct)(int,int,void*),
    void* matchFunctArguments) {
  wavefront_aligner_set_match_funct(wfAligner,matchFunct,matchFunctArguments);
}
void WFAligner::setMaxAlignmentScore(
    const int maxAlignmentScore) {
  wavefront_aligner_set_max_alignment_score(
      wfAligner,maxAlignmentScore);
}
void WFAligner::setMaxMemory(
    const uint64_t maxMemoryCompact,
    const uint64_t maxMemoryResident,
    const uint64_t maxMemoryAbort) {
  wavefront_aligner_set_max_memory(wfAligner,
      maxMemoryCompact,maxMemoryResident,maxMemoryAbort);
}
/*
 * Accessors
 */
int WFAligner::getAlignmentScore() {
  return wfAligner->cigar.score;
}
void WFAligner::getAlignmentCigar(
    char** const buffer,
    int* bufferLength) {
 *buffer = wfAligner->cigar.operations + wfAligner->cigar.begin_offset;
 *bufferLength = wfAligner->cigar.end_offset - wfAligner->cigar.begin_offset;
}
std::string WFAligner::getAlignmentCigar() {
  // Fetch CIGAR
  char* buffer;
  int bufferLength;
  getAlignmentCigar(&buffer,&bufferLength);
  // Create string and return
  return std::string(buffer,bufferLength);
}
/*
 * Align End-to-end
 */
int WFAligner::alignEnd2End(
    const int patternLength,
    const int textLength) {
  // Configure
  wavefront_aligner_set_alignment_end_to_end(wfAligner);
  // Align (using custom matching function)
  return wavefront_align(wfAligner,
      NULL,patternLength,
      NULL,textLength);
}
int WFAligner::alignEnd2End(
    const char* const pattern,
    const int patternLength,
    const char* const text,
    const int textLength) {
  // Configure
  wavefront_aligner_set_alignment_end_to_end(wfAligner);
  // Align
  return wavefront_align(wfAligner,
      pattern,patternLength,
      text,textLength);
}
int WFAligner::alignEnd2End(
    std::string& pattern,
    std::string& text) {
  // Delegate
  return alignEnd2End(
      pattern.c_str(),pattern.length(),
      text.c_str(),text.length());
}
/*
 * Align Ends-free
 */
// Align Ends-free
int WFAligner::alignEndsFree(
    const int patternLength,
    const int patternBeginFree,
    const int patternEndFree,
    const int textLength,
    const int textBeginFree,
    const int textEndFree) {
  // Configure
  wavefront_aligner_set_alignment_free_ends(wfAligner,
      patternBeginFree,patternEndFree,
      textBeginFree,textEndFree);
  // Align (using custom matching function)
  return wavefront_align(wfAligner,
      NULL,patternLength,
      NULL,textLength);
}
int WFAligner::alignEndsFree(
    const char* const pattern,
    const int patternLength,
    const int patternBeginFree,
    const int patternEndFree,
    const char* const text,
    const int textLength,
    const int textBeginFree,
    const int textEndFree) {
  // Configure
  wavefront_aligner_set_alignment_free_ends(wfAligner,
      patternBeginFree,patternEndFree,
      textBeginFree,textEndFree);
  // Align
  return wavefront_align(wfAligner,
      pattern,patternLength,
      text,textLength);
}
int WFAligner::alignEndsFree(
    std::string& pattern,
    const int patternBeginFree,
    const int patternEndFree,
    std::string& text,
    const int textBeginFree,
    const int textEndFree) {
  // Delegate
  return alignEnd2End(
      pattern.c_str(),pattern.length(),
      text.c_str(),text.length());
}
/*
 * Display & errors
 */
char* WFAligner::strError(
    const int wfErrorCode) {
  return wavefront_align_strerror(wfErrorCode);
}
/*
 * Gap-Affine Aligner
 */
WFAlignerGapAffine::WFAlignerGapAffine(
    const int mismatch,
    const int gapOpening,
    const int gapExtension,
    const bool onlyScore,
    const MemoryModel memoryModel) :
        WFAligner(onlyScore,memoryModel) {
  attributes.distance_metric = gap_affine;
  attributes.affine_penalties.match = 0;
  attributes.affine_penalties.mismatch = mismatch;
  attributes.affine_penalties.gap_opening = gapOpening;
  attributes.affine_penalties.gap_extension = gapExtension;
  wfAligner = wavefront_aligner_new(&attributes);
}
/*
 * Gap-Affine Dual-Cost (2 pieces) Aligner
 */
WFAlignerGapAffine2Pieces::WFAlignerGapAffine2Pieces(
    const int mismatch,
    const int gapOpening1,
    const int gapExtension1,
    const int gapOpening2,
    const int gapExtension2,
    const bool onlyScore,
    const MemoryModel memoryModel) :
        WFAligner(onlyScore,memoryModel) {
  attributes.distance_metric = gap_affine_2p;
  attributes.affine2p_penalties.match = 0;
  attributes.affine2p_penalties.mismatch = mismatch;
  attributes.affine2p_penalties.gap_opening1 = gapOpening1;
  attributes.affine2p_penalties.gap_extension1 = gapExtension1;
  attributes.affine2p_penalties.gap_opening2 = gapOpening2;
  attributes.affine2p_penalties.gap_extension2 = gapExtension2;
  wfAligner = wavefront_aligner_new(&attributes);
}



