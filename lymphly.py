__version__ = "1.0.0"

import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import warnings
from liftover import get_lifter
from tqdm import tqdm

config = {
    'feature_levels': ['Core', 'Extended'],
    'abc_subtypes_order': ['MCD', 'BN2', 'JS3', 'N1'],
    'gcb_subtypes_order': ['EZB', 'JS6'],
    'main_abc_subtypes': ['MCD', 'BN2'],
    'main_gcb_subtypes': ['EZB'],
    'dependent_statuses': ['MYC+', 'TP53+', 'A+']
}

cna_groups = {
    '10q23': ['FAS', 'PTEN'],
    '18q21': ['BCL2', 'MALT1'],
    '9p24': ['CD274', 'JAK2', 'PDCD1LG2'],
    '20q11': ['ASXL1', 'DNMT3B']
}

variant_classification_groups = {
    'MISSENSE': ['Missense_Mutation', 'Intron', "3'UTR", "5'UTR"],
    'NONSENSE': ['Nonsense_Mutation'],
    'FRAME_SHIFT': ['Frame_Shift_Ins', 'Frame_Shift_Del', 'Splice_Site']
}

palette_pure = {
    'EZB': '#1e90ff',
    'JS6': '#2ca02c',
    'JS3': '#9e0142',
    'BN2': '#ff7f0e',
    'MCD': '#d62728',
    'N1': '#b8860b',
    'MYC+': '#9467bd',
    'TP53+': '#7f7f7f',
    'A+': '#ff4500',
    'Other': '#c7c7c7'
}

def convert_hg19_to_hg38(input_maf: pd.DataFrame) -> pd.DataFrame:
    """
    Converts mutation positions from hg19 to hg38.
    Returns a new DataFrame with updated positions (hg38).
    """
    maf_conv = input_maf.copy()
    lifter_ = get_lifter('hg19', 'hg38')
    
    start_list = []
    end_list = []
    
    for i in tqdm(maf_conv.index, desc = '{Converting mutation positions from HG19 to HG38}'):
        chrom = maf_conv.loc[i, 'Chromosome']
        start = int(maf_conv.loc[i, 'Start_Position'])
        end = int(maf_conv.loc[i, 'End_Position'])

        # Convert start position
        lifted_start = lifter_[chrom][start]
        if lifted_start:
            start_list.append(int(lifted_start[0][1]))
        else:
            start_list.append(None)

        # Convert end position
        lifted_end = lifter_[chrom][end]
        if lifted_end:
            end_list.append(int(lifted_end[0][1]))
        else:
            end_list.append(None)

    # Update DataFrame
    maf_conv['Start_Position'] = start_list
    maf_conv['End_Position'] = end_list
    
    # Drop rows where either Start or End could not be converted
    maf_conv.dropna(subset=['Start_Position', 'End_Position'], inplace=True)
    
    maf_conv['Start_Position'] = maf_conv['Start_Position'].astype(int)
    maf_conv['End_Position'] = maf_conv['End_Position'].astype(int)
    
    if 'NCBI_Build' in maf_conv.columns:
        maf_conv['NCBI_Build'] = 'GRCh38'
    
    return maf_conv

def format_maf(input_maf: pd.DataFrame,
              variant_classifications: dict) -> pd.DataFrame:
    """
    Standardizes the MAF (Mutation Annotation Format) DataFrame by ensuring that necessary columns are present and correctly formatted.

    - Ensures that a valid column for sample IDs exists, either using 'Sample' or 'Tumor_Sample_Barcode'.
    - Computes and adds an 'End_Position' column if missing, based on the 'Start_Position' and 'Reference_Allele' columns.
    - Ensures the appropriate column for the alternative allele exists ('Tumor_Seq_Allele2' or 'Tumor_Seq_Allele1').
    - Formats the 'Chromosome' column to ensure it contains chromosome identifiers prefixed with 'chr'.
    - Checks for the presence of the 'Variant_Classification' column and updates 'Missense_Mutation' entries to 'Nonsense_Mutation' 
      if they have the mutation effect 'TRUNC'.

    Parameters:
    -----------
    input_maf : pandas.DataFrame
        A DataFrame containing mutation annotation data. It should include columns such as 'Sample', 'Chromosome', 
        'Start_Position', 'Reference_Allele', and others.

    Returns:
    --------
    pandas.DataFrame
        The modified DataFrame with standardized columns, including:
        - 'Sample' column for sample IDs (renamed from 'Tumor_Sample_Barcode' if necessary).
        - Computed 'End_Position' column if missing.
        - The 'Tumor_Seq_Allele2' column for the alternative allele.
        - Formatted 'Chromosome' column with 'chr' prefix.
        - Updates to the 'Variant_Classification' column where necessary.

    Raises:
    -------
    KeyError:
        If any required columns are missing, such as:
        - Neither 'Sample' nor 'Tumor_Sample_Barcode' are present.
        - 'Start_Position' or 'Reference_Allele' are missing when computing 'End_Position'.
        - 'Tumor_Seq_Allele2' or 'Tumor_Seq_Allele1' are missing for the alternative allele.
        - The 'Variant_Classification' column is missing.
    """
    maf = input_maf.copy()
    
    # Rename columns for consistency
    maf = maf.rename(columns={
        'Start_position': 'Start_Position',
        'End_position': 'End_Position',
    })
    
    # Check and rename 'Sample' column if needed
    if 'Sample' in maf.columns:
        print("'Sample' column in the MAF file will be used for sample IDs.")
    elif 'Tumor_Sample_Barcode' in maf.columns:
        print("There is no 'Sample' column in the MAF file, so the 'Tumor_Sample_Barcode' column will be used for sample IDs.")
        maf = maf.rename(columns={'Tumor_Sample_Barcode': 'Sample'})
    else:
        raise KeyError("There are no 'Sample' or 'Tumor_Sample_Barcode' columns in the MAF file. Please rename a column manually.")
    
    # Compute 'End_Position' if missing
    if 'End_Position' not in maf.columns:
        if 'Start_Position' in maf.columns and 'Reference_Allele' in maf.columns:
            print("There is no 'End_Position' column in the MAF file, so it will be added automatically.")
            maf['End_Position'] = maf['Start_Position'] + maf['Reference_Allele'].str.len() - 1
        else:
            raise KeyError("Cannot compute 'End_Position' because 'Start_Position' or 'Reference_Allele' is missing.")
    
    # Check and rename the alternative allele column
    if 'Tumor_Seq_Allele2' in maf.columns:
        print("'Tumor_Seq_Allele2' column will be used for alternative allele.")
    elif 'Tumor_Seq_Allele1' in maf.columns:
        print("'Tumor_Seq_Allele1' column will be used for alternative allele.")
        maf = maf.rename(columns={'Tumor_Seq_Allele1': 'Tumor_Seq_Allele2'})
    else:
        raise KeyError("There is no 'Tumor_Seq_Allele2' or 'Tumor_Seq_Allele1' column in the MAF file.")
    
    # Format the 'Chromosome' column to include 'chr'
    maf['Chromosome'] = maf['Chromosome'].astype(str).replace('^.*?(\\d+).*?$', 'chr\\1', regex = True)
    
    # Update 'Variant_Classification' for specific cases
    if 'Variant_Classification' in maf.columns:
        if 'Nonsense_Mutation' not in maf['Variant_Classification'].values:
            print("There are no 'Nonsense_Mutation' in 'Variant_Classification' column, so 'Missense_Mutation' with 'Mutation_effect'=='TRUNC' will be updated.")
            maf.loc[
                (maf['Variant_Classification'] == 'Missense_Mutation') & (maf['Mutation_effect'] == 'TRUNC'),
                'Variant_Classification'
            ] = 'Nonsense_Mutation'
    else:
        raise KeyError("The column 'Variant_Classification' is missing in the MAF file.")
    
    return maf

def check_cna_arm(cna_arm: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """
    Validates and filters a CNA arm matrix to ensure compatibility with annotation and feature tables.

    Parameters
    ----------
    cna_arm : pd.DataFrame
        DataFrame of arm-level CNA calls with chromosome arms as rows and samples as columns. 
        Expected values are integers in [0, 1, 2, -1, -2], representing deviations from sample ploidy.
    ann : pd.DataFrame
        Annotation table where the index represents sample identifiers.

    Returns
    -------
    pd.DataFrame
        Filtered version of `cna_arm` that:
        - Contains only valid chromosome arms `{1–24}{p/q}` (with 23 = X, 24 = Y). Raises an error if any invalid names are found.
        - Excludes arms (`23p`, `23q`, `24p`, `24q`, `21p`, `22p`, `13p`, `14p`, `15p`).
        - Contains only the remaining arms and all columns (samples) from the input matrix.
        - Raises a ValueError if unexpected CNA values are found.

    Raises
    ------
    ValueError
        - If `cna_arm` contains values outside the allowed CNA set `{0, 1, 2, -1, -2}`.
        - If any row index (chromosome arm) does not match the required `{1–24}{p/q}` format.

    Warnings
    --------
    - If samples in `ann` are missing from `cna_arm`, a warning is issued.
    - If extra samples are found in `cna_arm` but not in `ann`, they are dropped with an info message.
    """
    # Allowed CNA values
    allowed_values = [0, 1, 2, -1, -2]
    unique_vals = pd.unique(cna_arm.values.ravel()).tolist()
    unexpected_vals = [v for v in unique_vals if v not in allowed_values]
    if unexpected_vals:
        raise ValueError(f'Unexpected CNA values found in cna_arm: {unexpected_vals}')

    # Validate chromosome arm names
    valid_arms = [f"{i}{arm}" for i in range(1,25) for arm in ('p','q')]
    invalid_arms = [a for a in cna_arm.index if a not in valid_arms]
    if invalid_arms:
        raise ValueError(f'Invalid chromosome arms found: {invalid_arms}')

    # Exclude arms
    arms_to_exclude = ['23p','23q','24p','24q','21p','22p','13p','14p','15p']
    arms_to_keep = cna_arm.index.difference(arms_to_exclude)
    cna_arm_filtered = cna_arm.loc[arms_to_keep]

    # Filter columns based on annotation
    cna_cols = set(cna_arm_filtered.columns)
    ann_idx = set(ann.index)

    # Drop extra columns not in annotation
    extra_cols = cna_cols - ann_idx
    if extra_cols:
        print(f'Info: Columns in cna_arm not found among the annotation samples, these will be dropped: {sorted(extra_cols)}')
        cna_arm_filtered = cna_arm_filtered.drop(columns=extra_cols)

    # Warn about missing samples in cna_arm
    missing_cols = ann_idx - cna_cols
    if missing_cols:
        warnings.warn(f'Samples present in annotation but missing in cna_arm: {sorted(missing_cols)}')

    return cna_arm_filtered

def check_cna_gene(cna_gene: pd.DataFrame, ann: pd.DataFrame, feature_table: pd.DataFrame) -> pd.DataFrame:
    """
    Validates and filters a CNA gene matrix to ensure compatibility with annotation and feature tables.

    Parameters
    ----------
    cna_gene : pd.DataFrame
        DataFrame of CNA calls with genes as rows and samples as columns. 
        Expected values are in the [0, 1, 2, -1, -2].
    ann : pd.DataFrame
        Annotation table where the index represents sample identifiers.
    feature_table : pd.DataFrame
        Table listing features used for classification, including rows where Type == 'CNA_GENE'.

    Returns
    -------
    pd.DataFrame
        Filtered version of `cna_gene` that:
        - Contains only samples present in `ann`.
        - Emits a warning if expected genes (from `feature_table`) are missing.
        - Raises a ValueError if unexpected CNA values are found.

    Raises
    ------
    ValueError
        If `cna_gene` contains values outside of the allowed CNA set {0, 1, 2, -1, -2}.

    Warnings
    --------
    - If samples in `ann` are missing from `cna_gene`, a warning is issued.
    - If genes listed in the feature_table as 'CNA_GENE' are missing from `cna_gene`, a warning is issued.
    - If extra samples are found in `cna_gene` but not in `ann`, they are dropped with an info message.
    """
    allowed_values = [0, 1, 2, -1, -2]
    
    unique_vals = pd.unique(cna_gene.values.ravel()).tolist()
    unexpected_vals = [v for v in unique_vals if v not in allowed_values]
    if unexpected_vals:
        raise ValueError(f'Unexpected CNA values found in cna_gene: {unexpected_vals}')

    cna_cols = set(cna_gene.columns)
    ann_idx = set(ann.index)

    extra_cols = cna_cols - ann_idx
    if extra_cols:
        print(f'Info: Columns in cna_gene not found among the annotation samples, these will be dropped: {sorted(extra_cols)}')
        cna_gene_filtered = cna_gene.drop(columns=extra_cols)
    else:
        cna_gene_filtered = cna_gene.copy()

    missing_cols = ann_idx - cna_cols
    if missing_cols:
        warnings.warn(f'Columns in ann.index missing from cna_gene.columns: {sorted(missing_cols)}')

    cna_genes_in_table = set(feature_table.loc[feature_table.Type == 'CNA_GENE', 'Gene'].unique())
    cna_genes_in_cna_gene = set(cna_gene_filtered.index)
    missing_genes = cna_genes_in_table - cna_genes_in_cna_gene
    if missing_genes:
        warnings.warn(f'Genes present in feature_table CNA_GENE list but missing in cna_gene: {sorted(missing_genes)}')

    return cna_gene_filtered

def find_features(feature_attrs, variant_classifications, maf=None, cna_gene=None, cna_arm=None, sample_annot=None, features_to_samples={}, samples_seen_in_cna_group=None):
    """
    Identifies features in a mutation annotation file (MAF) based on a specific gene, its associated feature type, 
    and other relevant genomic data, and returns a list of features and associated samples.

    This function searches for different feature types:
    - Translocations (features marked as 'TRANSLOC').
    - Copy number alteration (CNA) events involving a specific gene (features marked as 'CNA_GENE').
    - Hotspots or mutations (missense, nonsense, frameshifts) within a specified genomic range (features marked as 'HOTSPOT').
    
    The function finds relevant mutations based on the gene, its chromosome, and genomic positions, and 
    associates them with the provided sample list. It ensures that no sample is counted more than once 
    for the same feature.

    Parameters:
    -----------
    maf : pandas.DataFrame
        The mutation annotation file (MAF) DataFrame containing mutation details, including gene names, 
        chromosome information, mutation types, and sample IDs.
        
    cna_gene : pandas.DataFrame
        A DataFrame containing copy number alterations (CNA) data, used to identify specific CNA events 
        for the given gene.
        
    feature_attrs : pandas.Series
        pd.Series containing information about the gene, its type (e.g., 'TRANSLOC', 'CNA_GENE', 'HOTSPOT'), 
        and its genomic coordinates.
        
    sample_annot : pandas.DataFrame
        A DataFrame containing annotation information for the samples, including whether certain conditions (e.g., 
        'True' for presence) are met for a specific gene and sample.
        
    features_to_samples: dict
        A dictionary tracking previously encountered features to ensure that samples are not double-counted 
        for the same feature.

    variant_classifications: dict
        A dictionary mapping variant classifications to feature groups
        
    Returns:
    --------
    list : [str, list, dict]
        A list containing:
        - A string representing the feature type (e.g., 'Gene_Transloc', 'Gene_HOTSPOT', etc.).
        - A list of sample IDs associated with the feature.
        - The updated `features_to_samples` dictionary with new sample information for the feature.
    """
    

    # Extract gene and feature type from feature_attrs
    gene, feat_type = feature_attrs[['Gene', 'Type']]
    
    # Handle translocations
    if feat_type == 'TRANSLOC':
        samples_with_feature = list(sample_annot[sample_annot[gene].isin(['True', 'TRUE', 'true', True])].index)
        feature_id = f'{gene}_{feat_type}'

    # Handle aneuploidy
    elif feat_type == 'ANEUPLOIDY':
        aneuploidy_treshold = int(feature_attrs['CNA_type'])
        nonzero_counts = (cna_arm != 0).sum(axis=0)
        samples_with_feature = nonzero_counts[nonzero_counts > aneuploidy_treshold].index.tolist()
        feature_id = f'{feat_type}'
        
    # Handle CNA gene events
    elif feat_type == 'CNA_GENE':
        if gene in cna_gene.index:
            if isinstance(feature_attrs['CNA_type'], str):
                cna_type = list(map(int, feature_attrs['CNA_type'].split(';')))
            else:
                cna_type = [int(feature_attrs['CNA_type'])]
            samples_with_feature = cna_gene.loc[:, cna_gene.loc[gene].isin(cna_type)].columns.tolist()

            group_name = None
            for gname, genes in cna_groups.items():
                if gene in genes:
                    group_name = gname
                    break

            if group_name is not None:
                new_samples = [s for s in samples_with_feature if s not in samples_seen_in_cna_group[group_name]]
                samples_seen_in_cna_group[group_name].update(new_samples)
                samples_with_feature = new_samples
                feature_id = f'{group_name}_{feat_type}'
            else:
                feature_id = f'{gene}_{feat_type}'
        else:
            samples_with_feature = []
            feature_id = f'{gene}_{feat_type}'
        
    # Handle hotspots or mutations (missense, nonsense, frameshifts)
    else:
        gene_chr_bool = (
            (maf['Hugo_Symbol'] == gene) &
            (maf['Chromosome'] == f"chr{int(feature_attrs['Chromosome'])}")
        )
        feat_start, feat_end = feature_attrs[['GRCh38_start', 'GRCh38_end']].astype(int)
        
        # Create search boolean mask based on feature type
        if feat_type == 'HOTSPOT':
            matching_maf_rows = (
                gene_chr_bool &
                (maf['Start_Position'] == feat_start) &
                (maf['End_Position'] == feat_end) &
                (maf['Reference_Allele'] == feature_attrs['Ref']) &
                (maf['Tumor_Seq_Allele2'] == feature_attrs['Alt'])
            )
            feature_id = '_'.join([gene, feature_attrs['Annotation']])
        else:
            matching_maf_rows = (
                gene_chr_bool &
                (maf['Start_Position'] > feat_start) &
                (maf['End_Position'] < feat_end) &
                (maf['Variant_Classification'].isin(variant_classifications.get(feat_type)))
            )
            feature_id = f'{gene}_{feat_start}-{feat_end}_{feat_type}'
        samples_with_feature = list()
        
        # Track samples and prevent duplicate entries for the same feature
        for _, mafrow in maf[matching_maf_rows].iterrows():
            feature_key_from_maf = '_'.join(mafrow[['Hugo_Symbol',
                                   'Chromosome', 
                                   'Start_Position',
                                   'End_Position',
                                   'Reference_Allele',
                                   'Tumor_Seq_Allele2']].astype(str)
                          )
            if feature_key_from_maf in features_to_samples:
                if mafrow['Sample'] not in features_to_samples[feature_key_from_maf]:
                    samples_with_feature.append(mafrow['Sample'])
                    features_to_samples[feature_key_from_maf].append(mafrow['Sample'])
            else:
                features_to_samples[feature_key_from_maf] = [mafrow['Sample']]
                samples_with_feature.append(mafrow['Sample'])
        
        # Remove duplicate samples
        samples_with_feature = list(set(samples_with_feature))
    
    return feature_id, samples_with_feature, features_to_samples

def finalize_subtype(lymphly_table, mutated_genes, subtypes, use_statuses, **kwargs):
    """
    Finalizes the subtype classification based on mutation data and genetic statuses.

    This function updates the 'Lymphly' column in the input 'lymphly_table' DataFrame based on various criteria:
    - Dominant core subtypes are prioritized based on unique gene mutations.
    - Extended subtypes are considered when core subtypes are absent.
    - The function supports the inclusion of genetic statuses if the 'use_statuses' flag is True.
    - It handles different levels of subtypes (core, extended) and includes specific logic for signature-based classifications.

    Parameters:
    -----------
    lymphly_table : pd.DataFrame
        A DataFrame with samples and their current subtype classification (if available) in the 'Lymphly' column.
    
    mutated_genes : dict
        A dictionary containing mutation data for each sample at different levels (Core, Extended).
        Format: {'Core': {sample: [list of mutated genes]}, 'Extended': {sample: [list of mutated genes]}}
    
    subtypes : dict
        A dictionary containing subtypes for each sample. Format: {level: {sample: [list of subtypes]}}
    
    use_statuses : bool
        A flag to determine whether to include genetic statuses in the final subtype classification.

    Returns:
    --------
    pd.DataFrame
        The updated 'lymphly_table' with the finalized subtype classifications in the 'Lymphly' column.
    
    Notes:
    ------
    The function first attempts to classify subtypes based on core and extended gene mutations.
    If no classification can be made based on genes, it will check for the presence of statuses.
    If no subtypes or statuses are found, the sample is classified as 'Other'.
    """
    abc_subtypes_order = kwargs.get('abc_subtypes_order')
    gcb_subtypes_order = kwargs.get('gcb_subtypes_order')
    main_abc_subtypes = kwargs.get('main_abc_subtypes')
    main_gcb_subtypes = kwargs.get('main_gcb_subtypes')
    dependent_statuses = kwargs.get('dependent_statuses')

    lymphly_table_subtype = lymphly_table.copy()
    cols = ['Lymphly'] + subtypes + ['Decision_Features', 'Algorithm_Step']

    for col in cols[::-1]:
        if col in subtypes:
            lymphly_table_subtype.insert(0, col, 'No')
        else:
            lymphly_table_subtype.insert(0, col, None)
    
    mutated_unique_genes_subtypes = {
        level: {
            sample: [
                subtype for gene in set(genes)
                if (subtype := gene.split('_', 1)[-1]) and '+' not in subtype
            ]
            for sample, genes in samples.items()
        }
        for level, samples in mutated_genes.items()
    }
    mutated_genes_subtypes = {
        level: {
            sample: [
                subtype for gene in genes
                if (subtype := gene.split('_', 1)[-1]) and '+' not in subtype
            ]
            for sample, genes in samples.items()
        }
        for level, samples in mutated_genes.items()
    }
    statuse_subtypes = {
        level: {
            sample: [
                subtype for gene in set(genes)
                if (subtype := gene.split('_', 1)[-1]) and '+' in subtype
            ]
            for sample, genes in samples.items()
        }
        for level, samples in mutated_genes.items()
    }

    ### Zero features to other
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_genes_subtypes['Core'][sample]
        extended_subtypes = mutated_genes_subtypes['Extended'][sample]
        core_statuses = statuse_subtypes['Core'][sample]
        extended_statuses = statuse_subtypes['Extended'][sample]
        statuses = core_statuses + extended_statuses
        all_subtypes = core_subtypes + extended_subtypes + statuses
        if not all_subtypes:
            lymphly_table_subtype.loc[sample, 'Lymphly'] = 'Other'
            lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'No features'

    ### Dominant core subtypes    
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_unique_genes_subtypes['Core'][sample]
        if core_subtypes:
            gene_to_subtype_count = pd.Series(core_subtypes).value_counts()
            if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                top_subtype = gene_to_subtype_count.idxmax()
                lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                decision_core = lymphly_table_subtype.loc[sample,'Core_mol_findings'][top_subtype]
                lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_core
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Core genes'
                lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'

    ### Dominant all subtypes (extended support core)
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_unique_genes_subtypes['Core'][sample]
        extended_subtypes = mutated_unique_genes_subtypes['Extended'][sample]
        all_subtypes = core_subtypes + extended_subtypes
        if core_subtypes:
            gene_to_subtype_count = pd.Series(all_subtypes).value_counts().loc[core_subtypes]
            if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                top_subtype = gene_to_subtype_count.idxmax()
                lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                decision_core = lymphly_table_subtype.loc[sample,'Core_mol_findings'][top_subtype]
                decision_extended = lymphly_table_subtype.loc[sample,'Extended_mol_findings'][top_subtype]
                lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_core+decision_extended
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Core+extended genes'
                lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'
                
    ### Dominant extended subtypes    
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_unique_genes_subtypes['Core'][sample]
        extended_subtypes = mutated_unique_genes_subtypes['Extended'][sample]
        if not core_subtypes:
            if extended_subtypes:
                gene_to_subtype_count = pd.Series(extended_subtypes).value_counts()
                if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                    top_subtype = gene_to_subtype_count.idxmax()
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                    decision_extended = lymphly_table_subtype.loc[sample,'Extended_mol_findings'][top_subtype]
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_extended
                    lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Extended genes'
                    lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'

    ### Dominant core subtypes (by mutations, not genes)
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_genes_subtypes['Core'][sample]
        if core_subtypes:
            gene_to_subtype_count = pd.Series(core_subtypes).value_counts()
            if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                top_subtype = gene_to_subtype_count.idxmax()
                lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                decision_core = lymphly_table_subtype.loc[sample,'Core_mol_findings'][top_subtype]
                lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_core
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Core mutations'
                lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'

    ### Dominant all subtypes (extended support core, by mutations, not genes)
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_genes_subtypes['Core'][sample]
        extended_subtypes = mutated_genes_subtypes['Extended'][sample]
        all_subtypes = core_subtypes + extended_subtypes
        if core_subtypes:
            gene_to_subtype_count = pd.Series(all_subtypes).value_counts().loc[core_subtypes]
            if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                top_subtype = gene_to_subtype_count.idxmax()
                lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                decision_core = lymphly_table_subtype.loc[sample,'Core_mol_findings'][top_subtype]
                decision_extended = lymphly_table_subtype.loc[sample,'Extended_mol_findings'][top_subtype]
                lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_core+decision_extended
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Core+extended mutations'
                lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'
    
    ### Dominant extended subtypes (by mutations, not genes)
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_genes_subtypes['Core'][sample]
        extended_subtypes = mutated_genes_subtypes['Extended'][sample]
        if not core_subtypes:
            if extended_subtypes:
                gene_to_subtype_count = pd.Series(extended_subtypes).value_counts()
                if (gene_to_subtype_count == gene_to_subtype_count.max()).sum() == 1:
                    top_subtype = gene_to_subtype_count.idxmax()
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = top_subtype
                    decision_extended = lymphly_table_subtype.loc[sample,'Extended_mol_findings'][top_subtype]
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = decision_extended
                    lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Extended mutations'
                    lymphly_table_subtype.loc[sample, top_subtype] = 'Yes'
    
    ### Presence among core features
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        core_subtypes = mutated_unique_genes_subtypes['Core'][sample]
        if core_subtypes:
            gene_to_subtype_count = pd.Series(core_subtypes).value_counts()
            max_count = gene_to_subtype_count.max()
            top_values = gene_to_subtype_count[gene_to_subtype_count == max_count].index
            core_subtypes = [x for x in core_subtypes if x in top_values]
            
            abc_present = [s for s in abc_subtypes_order if s in core_subtypes]
            gcb_present = [s for s in gcb_subtypes_order if s in core_subtypes]
            if abc_present and not gcb_present:
                selected_subtype = abc_present[0]
            elif gcb_present and not abc_present:
                selected_subtype = gcb_present[0]
            elif abc_present and gcb_present:
                abc_main = any(s in main_abc_subtypes for s in abc_present)
                gcb_main = any(s in main_gcb_subtypes for s in gcb_present)
                if abc_main and gcb_main:
                    selected_subtype = [gcb_present[0],abc_present[0]]
                elif abc_main and not gcb_main:
                    selected_subtype = abc_present[0]
                elif not abc_main and gcb_main:
                    selected_subtype = gcb_present[0]
                else:
                    selected_subtype = [gcb_present[0],abc_present[0]]
            if selected_subtype:
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Hierarchy core'
                lymphly_table_subtype.loc[sample, selected_subtype] = 'Yes'
                if isinstance(selected_subtype, list):
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = sum([lymphly_table_subtype.loc[sample, 'Core_mol_findings'][st] for st in selected_subtype], [])
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = '/'.join(sorted(selected_subtype))
                else:
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = lymphly_table_subtype.loc[sample,'Core_mol_findings'][selected_subtype]
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = selected_subtype

    ### Presence among extended features
    for sample in lymphly_table_subtype[lymphly_table_subtype.Lymphly.isna()].index:
        extended_subtypes = mutated_unique_genes_subtypes['Extended'][sample]
        if extended_subtypes:
            gene_to_subtype_count = pd.Series(extended_subtypes).value_counts()
            max_count = gene_to_subtype_count.max()
            top_values = gene_to_subtype_count[gene_to_subtype_count == max_count].index
            extended_subtypes = [x for x in extended_subtypes if x in top_values]
            
            abc_present = [s for s in abc_subtypes_order if s in extended_subtypes]
            gcb_present = [s for s in gcb_subtypes_order if s in extended_subtypes]
            if abc_present and not gcb_present:
                selected_subtype = abc_present[0]
            elif gcb_present and not abc_present:
                selected_subtype = gcb_present[0]
            elif abc_present and gcb_present:
                abc_main = any(s in main_abc_subtypes for s in abc_present)
                gcb_main = any(s in main_gcb_subtypes for s in gcb_present)
                if abc_main and gcb_main:
                    selected_subtype = [gcb_present[0],abc_present[0]]
                elif abc_main and not gcb_main:
                    selected_subtype = abc_present[0]
                elif not abc_main and gcb_main:
                    selected_subtype = gcb_present[0]
                else:
                    selected_subtype = [gcb_present[0],abc_present[0]]
            if selected_subtype:
                lymphly_table_subtype.loc[sample, 'Algorithm_Step'] = 'Hierarchy extended'
                lymphly_table_subtype.loc[sample, selected_subtype] = 'Yes'
                if isinstance(selected_subtype, list):
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = sum([lymphly_table_subtype.loc[sample, 'Extended_mol_findings'][st] for st in selected_subtype], [])
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = '/'.join(sorted(selected_subtype))
                else:
                    lymphly_table_subtype.at[sample, 'Decision_Features'] = lymphly_table_subtype.loc[sample,'Extended_mol_findings'][selected_subtype]
                    lymphly_table_subtype.loc[sample, 'Lymphly'] = selected_subtype

    ### Only statuses
    def process_only_statuses(row):
        core_subtypes = mutated_genes_subtypes['Core'][row.name]
        extended_subtypes = mutated_unique_genes_subtypes['Extended'][row.name]
        all_subtypes = core_subtypes + extended_subtypes
        statuses = statuse_subtypes['Extended'][row.name]
        if not all_subtypes and statuses:
            if not isinstance(row['Decision_Features'], list):
                row['Decision_Features'] = []
            row['Decision_Features'] = sum(
                [row['Extended_mol_findings'][st] for st in statuses], []
            )
            row['Lymphly'] = '/'.join(sorted(set(statuses)))
            row[statuses] = 'Yes'
            row['Algorithm_Step'] = 'Only statuses'
    
        return row

    lymphly_table_subtype = lymphly_table_subtype.apply(process_only_statuses, axis=1)

    ### Add statuses
    def update_decision_features(row):
        skip_subtypes = ['Other'] + list(dependent_statuses)
        if not any(s in skip_subtypes for s in row['Lymphly'].split('/')):
            row['Decision_Features'] = row['Decision_Features'].copy()
            statuses = statuse_subtypes['Extended'][row.name]
            if statuses:
                for status in statuses:
                    decision_extended = lymphly_table.loc[row.name, 'Extended_mol_findings'][status].copy()
                    row['Decision_Features'] += decision_extended
                    row[status] = 'Yes'
        return row
    
    if use_statuses:
        lymphly_table_subtype = lymphly_table_subtype.apply(update_decision_features, axis=1)
        
    def flatten_all_values(d):
        if isinstance(d, dict):
            combined = sum(d.values(), [])
            return '/'.join(combined)
        return d

    for col in ['Core_mol_findings', 'Extended_mol_findings']:
        lymphly_table_subtype[col] = lymphly_table_subtype[col].apply(flatten_all_values)    

    lymphly_table_subtype['Decision_Features'] = lymphly_table_subtype['Decision_Features'].apply(
        lambda x: '/'.join(list(set(x))) if isinstance(x, list) else x
    )
    
    return lymphly_table_subtype

def get_unique_subtype(feature_table):
    all_subtypes = set(feature_table.Subtype)
    no_plus = sorted([s for s in all_subtypes if '+' not in s])
    with_plus = sorted([s for s in all_subtypes if '+' in s])
    return no_plus + with_plus

def lymphly_classify(maf, cna_gene, cna_arm, annotation, feature_table, use_translocations, use_cna, use_statuses, **kwargs):
    """
    This function classifies samples based on mutation data and additional features like CNA,
    and generates a table with molecular findings and subtype assignments for each sample.

    Parameters:
    - maf (pd.DataFrame): Mutation data frame, including the 'Sample' column.
    - cna_gene (pd.DataFrame): CNA gene data for classified samples.
    - cna_arm (pd.DataFrame): CNA arm data for classified samples.
    - annotation (pd.DataFrame): Annotation table with predefined classifications.
    - feature_table (pd.DataFrame): Table of features used for classification, including mutation, CNA, etc.
    - use_translocations (bool): Whether to include translocations in the classification process.
    - use_cna (bool): Whether to include CNA data in the classification process.
    - use_statuses (bool): Whether to use genetic statuses for classification.

    Returns:
    - mutated_genes (dict): Dictionary of mutated genes for each sample and subtype level.
    - lymphly_table_subtype (pd.DataFrame): Data frame containing subtype assignments and molecular findings.
    - lymphly_statistics (dict): Dictionary of statistics for each classification level, including counts of features and subtypes.
    - subtypes (list): List of all subtypes used in the classification.
    """
    feature_levels = kwargs.get('feature_levels', [])
    abc_subtypes_order = kwargs.get('abc_subtypes_order')
    gcb_subtypes_order = kwargs.get('gcb_subtypes_order')
    main_abc_subtypes = kwargs.get('main_abc_subtypes')
    main_gcb_subtypes = kwargs.get('main_gcb_subtypes')
    dependent_statuses = kwargs.get('dependent_statuses')
    arms_to_exclude = kwargs.get('arms_to_exclude')

    # Get unique subtypes
    subtypes = get_unique_subtype(feature_table)
    features_to_samples = {}

    # Use only the samples from the annotation
    samples_to_classify = annotation.index.tolist()
    maf = maf[maf['Sample'].isin(samples_to_classify)]
    annotation = annotation.loc[samples_to_classify]

    # If translocations are not used, filter them out
    if not use_translocations:
        feature_table = feature_table[feature_table.Type != 'TRANSLOC']
        print('Translocations are not used in the classification.')
    else:
        print('Translocations are used in the classification.')
        translocs_not_in_annot = set(feature_table.loc[feature_table.Type == 'TRANSLOC', 'Gene']) - set(annotation.columns)
        if translocs_not_in_annot:
            warnings.warn(f"{translocs_not_in_annot} aren't in the annotation and will not be used for classification.")
            feature_table = feature_table.query('Type != "TRANSLOC" or Gene not in @translocs_not_in_annot')

    # If CNA is not used or cna_gene is None, filter out CNA features
    if not use_cna or cna_gene is None:
        feature_table = feature_table[~feature_table.Type.isin(['CNA_GENE'])]
        print('CNAs are not used in the classification.')
    else:
        print('CNAs are used in the classification.')
        cna_gene = cna_gene[[x for x in samples_to_classify if x in cna_gene.columns]]
        
    # If ARM is not used doesn't calculate aneuploidy status
    if cna_arm is None:
        feature_table = feature_table[~feature_table.Type.isin(['ANEUPLOIDY'])]
        print('Aneuploidy are not used in the classification.')
    else:
        print('Aneuploidy are used in the classification.')
        cna_arm = cna_arm[[x for x in samples_to_classify if x in cna_arm.columns]]

    # Create final table
    lymphly_table = maf.Sample.value_counts().to_frame('N_mutations')

    # Create list of mutated genes
    mutated_genes = {level: {sample: [] for sample in samples_to_classify} for level in feature_levels}

    # Process each level (Core/Extended)
    for level in feature_levels:

        samples_seen_in_cna_group = {group: set() for group in cna_groups}
        
        # Create table of features and subtypes
        feature_table_level = feature_table[feature_table.Level == level].copy()
        level_subtypes = get_unique_subtype(feature_table_level)
        lymphly_table_level = pd.DataFrame(0, index=samples_to_classify, columns=level_subtypes)
        lymphly_table_level[f'{level}_mol_findings'] = pd.Series(0, index=lymphly_table_level.index, dtype='object')

        # Iterative feature search
        for _, feature_attrs in feature_table_level.iterrows():
            feature_id, samples_with_feature, features_to_samples = find_features(feature_attrs, 
                                                                               variant_classification_groups,
                                                                               maf,
                                                                               cna_gene,
                                                                               cna_arm,
                                                                               annotation,
                                                                               features_to_samples,
                                                                               samples_seen_in_cna_group=samples_seen_in_cna_group
                                                                               )
            subtype = feature_attrs['Subtype']

            if samples_with_feature:
                sample_counts = pd.Series(samples_with_feature).value_counts()
                lymphly_table_level.loc[sample_counts.index, subtype] += sample_counts
                for sample, count in sample_counts.items():
                    mutated_genes[level][sample].append(f"{feature_attrs['Gene']}_{feature_attrs['Subtype']}")
                    current_dict = lymphly_table_level.at[sample, f'{level}_mol_findings']
                    if not isinstance(current_dict, dict):
                        current_dict = {}
                    subtype = feature_attrs['Subtype']
                    feature_entry = f"{feature_id}" if count == 1 else f"{count}*{feature_id}"
                    if subtype not in current_dict:
                        current_dict[subtype] = []
                    current_dict[subtype].append(feature_entry)
                    lymphly_table_level.at[sample, f'{level}_mol_findings'] = current_dict

        lymphly_table_level.rename(
            columns={col: f'{level}_{col}' for col in subtypes},
            inplace=True
        )
        lymphly_table_level[f'{level}_mol_findings'] = lymphly_table_level[f'{level}_mol_findings'].replace(0, '')
        lymphly_table = pd.concat([lymphly_table, lymphly_table_level], axis=1)

    # Finalize subtypes
    lymphly_table_subtype = finalize_subtype(lymphly_table, mutated_genes, subtypes, use_statuses, **kwargs)

    return lymphly_table_subtype

def pallete_preparation(lymphly_table, palette):
    unique_subtypes = set(lymphly_table.Lymphly)
    def mix_colors(colors):
        rgb_colors = [mcolors.to_rgb(c) for c in colors]
        avg_rgb = tuple(sum(x)/len(x) for x in zip(*rgb_colors))
        return mcolors.to_hex(avg_rgb)
    
    for subtype in unique_subtypes:
        parts = re.split(r'[-/]', subtype)
        if subtype not in palette_pure:
            if all(p in palette for p in parts):
                mixed = mix_colors([palette[p] for p in parts])
                palette[subtype] = mixed
            else:
                print(parts)
                palette[subtype] = '#bbbbbb'
    return palette

def subtypes_pieplot(lymphly_table, palette_all):
    subtype_series = pd.Series(lymphly_table.Lymphly)
    subtype_counts = subtype_series.value_counts().sort_index()
    colors = [palette_all.get(subtype, '#999999') for subtype in subtype_counts.index]
    
    def autopct_filter(pct):
        return f'{pct:.1f}%' if pct > 1.5 else ''
    
    def set_col_widths(table_obj, widths):
        for i, width in enumerate(widths):
            for key, cell in table_obj.get_celld().items():
                if key[1] == i:
                    cell.set_width(width)
    
    labels = [
        subtype if ('/' not in subtype and '+' not in subtype) or (subtype_counts[subtype] / subtype_counts.sum() * 100 > 1.5) else ''
        for subtype in subtype_counts.index
    ]
    
    table_data = pd.DataFrame({
        'Subtype': subtype_counts.index,
        'Count': subtype_counts.values,
        'Percentage': (subtype_counts / subtype_counts.sum() * 100).round(1)
    })
    all_row = pd.DataFrame({
        'Subtype': ['All'],
        'Count': [subtype_counts.sum()],
        'Percentage': [100.0]
    })
    table_data = pd.concat([table_data, all_row], ignore_index=True)
    
    if len(table_data) > 15:
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(36, 14))
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(25, 12))
        ax3 = None
    
    ax1.pie(
        subtype_counts,
        labels=labels,
        colors=colors,
        autopct=autopct_filter,
        startangle=90,
        counterclock=False,
        pctdistance=0.85,
        textprops={'fontsize': 22}
    )
    ax1.axis('equal')
    
    if ax3:
        half = len(table_data) // 2 + len(table_data) % 2
        part1 = table_data.iloc[:half]
        part2 = table_data.iloc[half:]
        
        ax2.axis('off')
        table2a = ax2.table(cellText=part1.values, colLabels=part1.columns, loc='center', cellLoc='center')
        table2a.set_fontsize(22)
        table2a.scale(1, 2)
        set_col_widths(table2a, [0.7, 0.2, 0.2])
    
        ax3.axis('off')
        table2b = ax3.table(cellText=part2.values, colLabels=part2.columns, loc='center', cellLoc='center')
        table2b.set_fontsize(22)
        table2b.scale(1, 2)
        set_col_widths(table2b, [0.7, 0.2, 0.2])
    else:
        ax2.axis('off')
        table = ax2.table(cellText=table_data.values, colLabels=table_data.columns, loc='center', cellLoc='center')
        table.set_fontsize(22)
        table.scale(1, 2)
        set_col_widths(table, [0.7, 0.2, 0.2])
    
    plt.tight_layout()
    plt.show()
