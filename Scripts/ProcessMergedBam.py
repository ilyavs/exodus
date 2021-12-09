import pysam
import pandas as pd
import numpy as np


def process_merged_bam(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as samfile:
        genome_sizes = pd.DataFrame(zip(samfile.references, samfile.lengths), columns=['ref_id', 'length']).set_index('ref_id').squeeze()
        data = []
        for a in samfile:
            if not a.is_unmapped:
                data.append({'read_name': a.query_name, 'is_read1': a.is_read1, 'read_length': a.query_length, 'alignment_length': a.query_alignment_length, 'is_supplementary': a.is_supplementary, 'ref_id': a.reference_name, 'alignment_score': a.get_tag('AS'), 'is_reverse': a.is_reverse, 'start': a.reference_start, 'end': a.reference_end, 'tlen': a.template_length, 'is_proper_pair': a.is_proper_pair})
    df = pd.DataFrame(data)
    df.set_index(['read_name', 'is_read1', 'ref_id'], inplace=True)
    df.drop(df.query('is_supplementary == True').index, inplace=True)
    df.drop('is_supplementary', axis=1, inplace=True)
    df = df[df.read_length == df.alignment_length] # filter clipped alignemnts, they won't get perfect score anyway
    dtypes = df.dtypes.to_dict()
    df = df.unstack('is_read1').dropna().stack('is_read1') # keeps only reads that both R1 and R2 mapped
    df = df.astype(dtypes) # restore dtypes
    orientation_bool = df.is_reverse.unstack('is_read1').sum(axis=1) == 1 # FR or RF orientation is True
    df = df.unstack('is_read1').loc[orientation_bool[orientation_bool].index].stack('is_read1') # Keep only reads in FR and RF orientation
    # Set orientation for each read. FR is when is_reverse = False and start <= end of mate. RF is when is_reverse = False and start > end of mate.
    df_reshaped = df.droplevel('is_read1').set_index('is_reverse', append=True).unstack('is_reverse')
    df_reshaped[('RF', False)] = False
    df_reshaped[('RF', True)] = False
    df_reshaped.loc[df_reshaped[('start', False)] > df_reshaped[('end', True)], ('RF', False)] = True
    df_reshaped.loc[df_reshaped[('start', False)] > df_reshaped[('end', True)], ('RF', True)] = True
    # Calculate corrected template length for RF read pairs
    max_proper_tlen = df[df.is_proper_pair].tlen.abs().max()
    df_reshaped_RF = df_reshaped[df_reshaped[('RF', False)]].reset_index()
    corrected_tlen = genome_sizes.reindex(df_reshaped_RF.ref_id).reset_index()['length'] - df_reshaped_RF[('start', False)] + df_reshaped_RF[('end'), True]
    df_reshaped_RF[('tlen', False)] = corrected_tlen
    df_reshaped_RF[('tlen', True)] = -1 * corrected_tlen
    df_reshaped_RF.set_index(df_reshaped.index.names, inplace=True)
    df_reshaped.update(df_reshaped_RF)
    # Filter for too large template length
    df_reshaped = df_reshaped[df_reshaped[('tlen', False)] <= max_proper_tlen]
    # Filter for too small template length
    df_reshaped = df_reshaped[df_reshaped[('tlen', False)] >= 0]
    df = df.unstack('is_read1').reindex(df_reshaped.index).stack('is_read1')
    # Filter same best ref
    read_lengths = df['read_length'].squeeze().unstack('ref_id').fillna(method='bfill', axis=1).fillna(method='ffill', axis=1)
    alignment_scores = df['alignment_score'].squeeze().unstack('ref_id').fillna(0)
    temp = (alignment_scores == pd.DataFrame(np.tile(alignment_scores.max(axis=1), (alignment_scores.shape[1], 1)).T, index=alignment_scores.index, columns=alignment_scores.columns)).unstack('is_read1').stack('ref_id').all(axis=1).unstack('ref_id').any(axis=1)
    alignment_scores = alignment_scores.unstack('is_read1').loc[temp[temp].index].stack('is_read1')
    read_lengths = read_lengths.unstack('is_read1').loc[temp[temp].index].stack('is_read1')
    return alignment_scores, read_lengths, genome_sizes


def process_reads_df(alignment_scores, read_lengths, genome_sizes):
    alignment_scores_reshaped = alignment_scores.unstack('is_read1')
    read_lengths_reshaped = read_lengths.unstack('is_read1')
    singles_index = (alignment_scores_reshaped.stack('ref_id') == read_lengths_reshaped.stack('ref_id')).all(axis=1).unstack('ref_id').sum(axis=1) == 1
    singles_scores = alignment_scores_reshaped[singles_index]
    singles_read_lengths = read_lengths_reshaped[singles_index]
    genome_counts_from_singles = (singles_scores == singles_read_lengths).stack('ref_id').all(axis=1).unstack('ref_id').sum() * 2
    bases = singles_read_lengths[(singles_scores == singles_read_lengths).stack('ref_id').all(axis=1).unstack('ref_id')].stack('is_read1').sum().astype(int)
    coverages = bases / genome_sizes
    ties_index = (alignment_scores_reshaped.stack('ref_id') == read_lengths_reshaped.stack('ref_id')).all(axis=1).unstack('ref_id').sum(axis=1) > 1
    if ties_index.sum() > 0:
        ties_scores = alignment_scores_reshaped[ties_index]
        ties_read_lengths = read_lengths_reshaped[ties_index]
        ties_bool_index = (ties_scores == ties_read_lengths).stack('ref_id').all(axis=1).unstack('ref_id')
        ties_coverages = ties_bool_index.replace(True, coverages).replace(False, 0)
        ties_coverages = ties_coverages[ties_coverages.sum(axis=1) > 0]
        ties_relative_abundances = ties_coverages.divide(ties_coverages.sum(axis=1), axis=0)
        # Fix ties for cases where there is only one genome with coverage
        selected_ties_relative_abundances = ties_relative_abundances[ties_relative_abundances.max(1) == 1]
        if selected_ties_relative_abundances.shape[0] > 0:
            selected_ties_read_lengths = ties_read_lengths.loc[ties_relative_abundances[ties_relative_abundances.max(1) == 1].index]
            ties_assignments = selected_ties_relative_abundances.apply(lambda x: x.index[x.astype('boolean')].to_frame().squeeze(), axis=1)
            genome_counts_from_ties = ties_assignments.value_counts()
            final_counts = pd.merge(genome_counts_from_singles.to_frame('counts_from_singles'),
                                    genome_counts_from_ties.to_frame('counts_from_ties'), left_index=True, right_index=True,
                                    how='left').fillna(0).sum(1)
            temp = ties_assignments.reset_index().rename(columns={0: 'ref_id'})
            temp['Value'] = 1
            ties_bases = (temp.set_index(['read_name', 'ref_id']).squeeze(axis=1).unstack('ref_id').reindex(columns=selected_ties_relative_abundances.columns).fillna(0) * selected_ties_read_lengths.stack('ref_id').sum(axis=1).unstack('ref_id')).sum()
            final_bases = bases + ties_bases
            final_coverages = final_bases / genome_sizes
            if ties_relative_abundances[ties_relative_abundances.max(1) < 1].shape[0] > 0:
                # Adjust the coverages
                ties_read_lengths = read_lengths_reshaped.loc[ties_relative_abundances[ties_relative_abundances.max(1) < 1].index]
                ties_bool_index = ties_bool_index.loc[ties_relative_abundances[ties_relative_abundances.max(1) < 1].index]
                ties_coverages = ties_bool_index.replace(True, final_coverages).replace(False, 0)
                ties_coverages = ties_coverages[ties_coverages.sum(axis=1) > 0]
                ties_relative_abundances = ties_coverages.divide(ties_coverages.sum(axis=1), axis=0)
                # Fix the rest of the ties
                rng = np.random.default_rng(42)
                ties_assignments = ties_relative_abundances.apply(lambda x: rng.choice(x.index, p=x), axis=1)
                genome_counts_from_ties = ties_assignments.value_counts()
                final_counts = pd.merge(final_counts.to_frame('counts_from_selected_ties'),
                                        genome_counts_from_ties.to_frame('counts_from_ties'), left_index=True, right_index=True,
                                        how='left').fillna(0).sum(1)
                temp = ties_assignments.reset_index().rename(columns={0: 'ref_id'})
                temp['Value'] = 1
                ties_bases = (temp.set_index(['read_name', 'ref_id']).squeeze(axis=1).unstack('ref_id').reindex(columns=ties_relative_abundances.columns).fillna(0) * ties_read_lengths.stack('ref_id').sum(axis=1).unstack('ref_id')).sum()
                final_bases = final_bases + ties_bases
        else:
            rng = np.random.default_rng(42)
            ties_assignments = ties_relative_abundances.apply(lambda x: rng.choice(x.index, p=x), axis=1)
            genome_counts_from_ties = ties_assignments.value_counts()
            final_counts = pd.merge(genome_counts_from_singles.to_frame('counts_from_singles'), genome_counts_from_ties.to_frame('counts_from_ties'), left_index=True, right_index=True, how='left').fillna(0).sum(1)
            temp = ties_assignments.reset_index().rename(columns={0: 'ref_id'})
            temp['Value'] = 1
            ties_bases = (temp.set_index(['read_name', 'ref_id']).squeeze(axis=1).unstack('ref_id').reindex(columns=ties_relative_abundances.columns).fillna(0) * ties_read_lengths.stack('ref_id').sum(axis=1).unstack('ref_id')).sum()
            final_bases = bases + ties_bases
    else:
        final_counts = genome_counts_from_singles
        final_bases = bases
    final_coverages = final_bases / genome_sizes
    final_relative_abundances = final_coverages/final_coverages.sum()
    return final_relative_abundances


# bam_path = '/data/workarea/ilyav/PipelineRuns/SampleOriented/ExodusCodeReview/20Nov829/exodus/1.align/20Nov829_bwa.merged.sorted.bam'
bam_path = snakemake.input[0]
relative_abundances_path = snakemake.output.relative_abundances_path

alignment_scores, read_lengths, genome_sizes = process_merged_bam(bam_path)
final_relative_abundances = process_reads_df(alignment_scores, read_lengths, genome_sizes)
final_relative_abundances = final_relative_abundances.to_frame().T
final_relative_abundances.insert(0, 'Sample', snakemake.wildcards.sample)
final_relative_abundances.to_csv(relative_abundances_path, index=False)
