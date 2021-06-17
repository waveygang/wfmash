// Configure the attributes of the wf-aligner
wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
attributes.distance_metric = gap_affine;
attributes.affine_penalties = affine_penalties;
// attributes.distance_metric = gap_affine2p;
// attributes.affine2p_penalties = affine2p_penalties;
attributes.reduction.reduction_strategy = wavefront_reduction_none; // wavefront_reduction_dynamic
attributes.reduction.min_wavefront_length = 10;
attributes.reduction.max_distance_threshold = 50;
attributes.alignment_scope = alignment_scope_alignment; // alignment_scope_score
attributes.low_memory = false;
wavefront_aligner_t* const wf_aligner =
    wavefront_aligner_new(strlen(pattern),strlen(text),&attributes);
// Clear & resize
wavefront_aligner_t* const wf_aligner = align_input->wf_aligner;
wavefront_aligner_clear__resize(
    wf_aligner,
    align_input->pattern_length,align_input->text_length);
// Align
wavefront_align(
    wf_aligner,
    pattern,strlen(pattern),text,strlen(text));
// CIGAR
cigar_print_pretty(
    stderr,
    pattern,strlen(pattern),text,strlen(text),
    &wf_aligner->cigar,mm_allocator);
fprintf(stderr,"SCORE: %d \n",
        cigar_score_gap_affine2p(&wf_aligner->cigar,&affine2p_penalties));
// Free
wavefront_aligner_delete(wf_aligner);
