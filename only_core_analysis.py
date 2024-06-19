# Date: 061924
import pandas as pd
import argparse
import os
import warnings
warnings.filterwarnings('ignore')


def get_oc_from_roary_gpa(file_path, cols_of_interest, cols_other):
    """
    Parse the gpa.csv file to extract only-core of group
    """
    # Read the CSV file
    df = pd.read_csv(file_path, sep=",")
    
    # Filter rows based on the presence in the specified columns and absence in the other columns
    filtered_df = df[
        df[cols_of_interest[0]].notna() & df[cols_of_interest[0]].astype(bool) & 
        df[cols_of_interest[1]].notna() & df[cols_of_interest[1]].astype(bool) & 
        df[cols_other].isna().all(axis=1)
    ]
    
    return filtered_df


def reorder_filtered_gpa(filtered_gpa_df, custom_cols):
    """
    Reorder and get only the columns of interest
    """
    reordered = filtered_gpa_df[custom_cols]
    return reordered


def best_hit_pgap(df, gpa_list, pident_cut, eval_cut):
    """
    Find the best hit for each query of Blastp, with Pident and e-value cutoff
    Note that gpa_list should be shared between Groups - "qseqid" is common
    """
    # Define the criteria for filtering
    criteria = (df['pident'] >= pident_cut) & (df['evalue'] < eval_cut)

    # Apply the criteria to filter the DataFrame
    filtered_df = df[criteria]

    # Sort by qseqid and bitscore to ensure the best hit is first for each qseqid
    sorted_df = filtered_df.sort_values(by=['qseqid', 'bitscore'], ascending=[True, False])

    # Select the best hit for each qseqid
    best_hits = sorted_df.drop_duplicates(subset='qseqid', keep='first')

    # Find terms that are missing
    best_hits_queries = list(best_hits['qseqid'])
    miss_list = []
    for term in gpa_list:
        if term not in best_hits_queries:
            miss_list.append(term)

    # Add missing terms to best_hits
    for miss in miss_list:
        blank_row = [miss] + [""] * 12
        best_hits.loc[-1] = blank_row
        best_hits.index = best_hits.index + 1
        best_hits = best_hits.sort_index()

    # Sort the best_hits based on gpa order
    best_hits['qseqid'] = pd.Categorical(best_hits['qseqid'], categories=gpa_list, ordered=True)
    best_hits = best_hits.sort_values(by='qseqid')

    return best_hits


def subset_meta(meta, lt_list):
    """
    Subset metadata with current LocusTag to match annotation and reorder them
    """
    lt_to_meta = meta[meta["LocusTag"].isin(lt_list)]
    cat_lt = []
    for l in lt_list:
        if l in cat_lt:
            pass
        else:
            cat_lt.append(l)
    lt_to_meta['LocusTag'] = pd.Categorical(lt_to_meta['LocusTag'], categories=cat_lt, ordered=True)
    lt_to_meta = lt_to_meta.sort_values('LocusTag')

    return lt_to_meta


def subset_meta2(meta, lt_list):
    """
    Subset metadata with current LocusTag to match annotation and reorder them,
    allowing for duplicates and blanks in lt_list
    """
    # Create an empty DataFrame to store the ordered results
    ordered_meta = pd.DataFrame()

    # Iterate through lt_list and extract corresponding rows from meta
    for lt in lt_list:
        if lt == "":
            # Append a blank row if lt is empty
            blank_row = pd.Series({col: "" for col in meta.columns})
            ordered_meta = pd.concat([ordered_meta, pd.DataFrame([blank_row])], ignore_index=True)
        else:
            matching_rows = meta[meta["LocusTag"] == lt]
            ordered_meta = pd.concat([ordered_meta, matching_rows], ignore_index=True)
    
    return ordered_meta


def best_hit_oc_nf(df, gpa_list, pident_cut, eval_cut):
    """
    Function to find the best hit for each query of Blastp
    Slightly modified version for Only-core-nf search
    """
    # Get all unique query ids in order
    query = list(df["qseqid"])
    unique_query = []
    for q in query: 
        if q not in unique_query: unique_query.append(q)

    # Define the criteria for filtering
    criteria = (df['pident'] >= pident_cut) & (df['evalue'] < eval_cut)

    # Apply the criteria to filter the DataFrame
    filtered_df = df[criteria]

    # Sort by qseqid and bitscore to ensure the best hit is first for each qseqid
    sorted_df = filtered_df.sort_values(by=['qseqid', 'bitscore'], ascending=[True, False])

    # Select the best hit for each qseqid
    best_hits = sorted_df.drop_duplicates(subset='qseqid', keep='first')

    # Check if Only-core-nf are present
    best_query = list(best_hits['qseqid'])
    for q in unique_query:
        if q not in best_query:
            # Insert a blank row to indiciate nf
            blank_row = [q] + [""] * 12
            best_hits.loc[-1] = blank_row
            best_hits.index = best_hits.index + 1
            best_hits = best_hits.sort_index()
    
    # Sort the dataframe
    best_hits['qseqid'] = pd.Categorical(best_hits['qseqid'], categories=gpa_list, ordered=True)
    best_hits = best_hits.sort_values(by='qseqid')

    # Reset index
    best_hits = best_hits.reset_index(drop=True)

    return best_hits


def main():
    # Argument parser
    parser = argparse.ArgumentParser(
        description="""
        Finding Only-cores and Only-core-NFs in 2 : 2 setting.\n
        This code assumes two groups with two strains each.\n
        Strain #1 and #2 are Group #1 strains. Strain #3 and #4 are Group #2 strains.\n
        +------------------- Requirement files before running -------------------+
        0. gene_presence_absence.csv | Roary result of merged groups (4 strains)
        1. Meta | LocusTag : Description mapping table of NCBI PGAP annotation
        2. Group info | This should match column name of gene_presence_absence.csv
        3. Blastp results of OCs to its PGAP annotation (X4 sets)
        4. Blastp results of OCs to other Group Strain (X4 sets)
        +------------------------------------------------------------------------+\n
        """,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w', '--workspace', required=True, help='Path to workspace path')
    parser.add_argument('-r', '--refspace', required=True, help='Path to reference (meta) space path')

    parser.add_argument('-g', '--gpa', required=True, help='Path to gene_presence_absence.csv')
    parser.add_argument('-m1', '--meta1', required=True, help='Path to meta.csv of Group #1 Strain #1')
    parser.add_argument('-m2', '--meta2', required=True, help='Path to meta.csv of Group #1 strain #2')
    parser.add_argument('-m3', '--meta3', required=True, help='Path to meta.csv of Group #2 strain #3')
    parser.add_argument('-m4', '--meta4', required=True, help='Path to meta.csv of Group #2 strain #4')

    parser.add_argument('-gi1', '--ginfo1', required=True, help='Name of Group #1', action='append')
    parser.add_argument('-gi2', '--ginfo2', required=True, help='Name of Group #2', action='append')

    parser.add_argument('-p1', '--pgap1', required=True, help='Path to Blastp of G1 OC-S1 to PGAP')
    parser.add_argument('-p2', '--pgap2', required=True, help='Path to Blastp of G1 OC-S2 to PGAP')
    parser.add_argument('-p3', '--pgap3', required=True, help='Path to Blastp of G2 OC-S3 to PGAP')
    parser.add_argument('-p4', '--pgap4', required=True, help='Path to Blastp of G2 OC-S4 to PGAP')

    parser.add_argument('--g1s3', required=True, help='Path to Blastp of G1 OC to S3')
    parser.add_argument('--g1s4', required=True, help='Path to Blastp of G1 OC to S4')
    parser.add_argument('--g2s1', required=True, help='Path to Blastp of G2 OC to S1')
    parser.add_argument('--g2s2', required=True, help='Path to Blastp of G2 OC to S2')

    parser.add_argument('--pident_cut_pgap', required=False, default=35, help='Pident cutoff for PGAP')
    parser.add_argument('--evalue_cut_pgap', required=False, default=0.001, help='Evalue cutoff for PGAP')
    parser.add_argument('--pident_cut_nf', required=False, default=35, help='Pident cutoff for NF')
    parser.add_argument('--evalue_cut_nf', required=False, default=0.001, help='Evalue cutoff for NF')

    parser.add_argument('-s', '--save_path', required=False, default='./', help='Path to save final data')
    parser.add_argument('--save_name1', required=False, default='group1_final.csv', help='G1 final name')
    parser.add_argument('--save_name2', required=False, default='group2_final.csvs', help='G2 final name')

    args = parser.parse_args() # Made by HY

    # Define argument variables
    workspace = args.workspace
    refspace = args.refspace

    gpa = args.gpa
    meta1 = args.meta1
    meta2 = args.meta2
    meta3 = args.meta3
    meta4 = args.meta4

    group_info1 = args.ginfo1
    group_info2 = args.ginfo2

    pgap1 = args.pgap1
    pgap2 = args.pgap2
    pgap3 = args.pgap3
    pgap4 = args.pgap4

    g1s3 = args.g1s3
    g1s4 = args.g1s4
    g2s1 = args.g2s1
    g2s2 = args.g2s2

    pident_cut_pgap = args.pident_cut_pgap
    evalue_cut_pgap = args.evalue_cut_pgap
    pident_cut_nf = args.pident_cut_nf
    evalue_cut_nf = args.evalue_cut_nf

    save_path = args.save_path
    save_name1 = args.save_name1
    save_name2 = args.save_name2

    # Run
    # -------------------------------------------------------------------------------------------------
    # 1. Define Only-core from gene_presence_absence.csv
    # -------------------------------------------------------------------------------------------------
    ## Path to gpa
    gpa_file_path = workspace + gpa

    ## Group #1 (GrowthSuppressed)
    g1_filtered_gpa = get_oc_from_roary_gpa(file_path=gpa_file_path,
                                            cols_of_interest=group_info1,
                                            cols_other=group_info2)
    print("Number of Group #1 only-core(s):", len(g1_filtered_gpa))

    ## Group #2 (NonGS)
    g2_filtered_gpa = get_oc_from_roary_gpa(file_path=gpa_file_path,
                                            cols_of_interest=group_info2,
                                            cols_other=group_info1)
    print("Number of Group #2 only-core(s):", len(g2_filtered_gpa))

    ## Group #1 | ID_1 : ID_2 : Gene : Non-unique Gene name : Prokka_Annotation
    g1_custom_cols = group_info1 + ["Gene", "Non-unique Gene name", "Annotation"]
    g1_filt_gpa_reordered = reorder_filtered_gpa(filtered_gpa_df=g1_filtered_gpa,
                                                 custom_cols=g1_custom_cols)
    
    ## Group #1 | Rename columns
    rename = ["ID " + "[" + x + "]" for x in group_info1]
    rename = rename + ["Gene", "Non-unique Gene name", "Annotation [By Prokka]"]
    g1_filt_gpa_reordered = g1_filt_gpa_reordered.set_axis(rename, axis=1)
    g1_filt_gpa_reordered = g1_filt_gpa_reordered.reset_index(drop=True)
    
    ## Group #2 | ID_1 : ID_2 : Gene : Non-unique Gene name : Prokka_Annotation
    g2_custom_cols = group_info2 + ["Gene", "Non-unique Gene name", "Annotation"]
    g2_filt_gpa_reordered = reorder_filtered_gpa(filtered_gpa_df=g2_filtered_gpa,
                                                 custom_cols=g2_custom_cols)
    
    ## Group #2 | Rename columns
    rename = ["ID " + "[" + x + "]" for x in group_info2]
    rename = rename + ["Gene", "Non-unique Gene name", "Annotation [By Prokka]"]
    g2_filt_gpa_reordered = g2_filt_gpa_reordered.set_axis(rename, axis=1)
    g2_filt_gpa_reordered = g2_filt_gpa_reordered.reset_index(drop=True)

    # -------------------------------------------------------------------------------------------------
    # 2. Prokka : PGAP Matching
    ## Prokka IDs are mapped to PGAP-annotated LocusTags
    ## LocusTags are then mapped to PGAP Gene annotation
    # -------------------------------------------------------------------------------------------------
    # Read PGAP meta files
    meta1 = pd.read_csv(refspace + meta1)
    meta2 = pd.read_csv(refspace + meta2)
    meta3 = pd.read_csv(refspace + meta3)
    meta4 = pd.read_csv(refspace + meta4)

    # Blastp Header
    blastp_header = ["qseqid","sseqid", "pident", "length", "mismatch", 
                     "gapopen", "qstart", "qend", "sstart", "send", 
                     "evalue", "bitscore", "qcovs"]
    
    # Group #1 | OC to PGAP matching Blastp results
    ## Strain #1 
    g1_s1_pgap = pd.read_csv(workspace + pgap1, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits for PGAP matching
    s1_gpa_list = g1_filtered_gpa[group_info1[0]]
    g1_s1_pgap_best = best_hit_pgap(df=g1_s1_pgap, gpa_list=s1_gpa_list,
                                    pident_cut=pident_cut_pgap, eval_cut=evalue_cut_pgap)

    ## LocusTag : PGAP Annotation
    s1_lt_to_annot = subset_meta2(meta=meta1, lt_list=list(g1_s1_pgap_best["sseqid"]))
    rename = ["LocusTag " + "[" + group_info1[0] + "]",
              "PGAP Gene annotation " + "[" + group_info1[0] + "]"]
    s1_lt_to_annot = s1_lt_to_annot.set_axis(rename, axis=1)
    s1_lt_to_annot = s1_lt_to_annot.reset_index(drop=True)

    ## Strain #2
    g1_s2_pgap = pd.read_csv(workspace + pgap2, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits for PGAP matching
    s2_gpa_list = g1_filtered_gpa[group_info1[1]]
    g1_s2_pgap_best = best_hit_pgap(df=g1_s2_pgap, gpa_list=s1_gpa_list,
                                    pident_cut=pident_cut_pgap, eval_cut=evalue_cut_pgap)

    ## LocusTag : PGAP Annotation
    s2_lt_to_annot = subset_meta2(meta=meta2, lt_list=list(g1_s2_pgap_best["sseqid"]))
    rename = ["LocusTag " + "[" + group_info1[1] + "]",
              "PGAP Gene annotation " + "[" + group_info1[1] + "]"]
    s2_lt_to_annot = s2_lt_to_annot.set_axis(rename, axis=1)
    s2_lt_to_annot = s2_lt_to_annot.reset_index(drop=True)

    # Group #1 | Merge LT : PGAP annotation dataframes
    g1_lt_to_annot = pd.concat([s1_lt_to_annot, s2_lt_to_annot], axis=1)

    # Group #1 | Merge Prokka & PGAP annotations
    g1_prokka_pgap = pd.concat([g1_filt_gpa_reordered, g1_lt_to_annot], axis=1)

    #

    # Group #2 | OC to PGAP matching Blastp results
    ## Strain #3
    g2_s3_pgap = pd.read_csv(workspace + pgap3, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits for PGAP matching
    s3_gpa_list = g2_filtered_gpa[group_info2[0]]
    g2_s3_pgap_best = best_hit_pgap(df=g2_s3_pgap, gpa_list=s3_gpa_list,
                                    pident_cut=pident_cut_pgap, eval_cut=evalue_cut_pgap)

    ## LocusTag : PGAP Annotation
    s3_lt_to_annot = subset_meta2(meta=meta3, lt_list=list(g2_s3_pgap_best["sseqid"]))
    rename = ["LocusTag " + "[" + group_info2[0] + "]",
              "PGAP Gene annotation " + "[" + group_info2[0] + "]"]
    s3_lt_to_annot = s3_lt_to_annot.set_axis(rename, axis=1)
    s3_lt_to_annot = s3_lt_to_annot.reset_index(drop=True)

    ## Strain #4
    g2_s4_pgap = pd.read_csv(workspace + pgap4, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits for PGAP matching
    s4_gpa_list = g2_filtered_gpa[group_info2[1]]
    g2_s4_pgap_best = best_hit_pgap(df=g2_s4_pgap, gpa_list=s3_gpa_list,
                                    pident_cut=pident_cut_pgap, eval_cut=evalue_cut_pgap)

    ## LocusTag : PGAP Annotation
    s4_lt_to_annot = subset_meta2(meta=meta4, lt_list=list(g2_s4_pgap_best["sseqid"]))
    rename = ["LocusTag " + "[" + group_info2[1] + "]",
              "PGAP Gene annotation " + "[" + group_info2[1] + "]"]
    s4_lt_to_annot = s4_lt_to_annot.set_axis(rename, axis=1)
    s4_lt_to_annot = s4_lt_to_annot.reset_index(drop=True)

    # Group #2 | Merge LT : PGAP annotation dataframes
    g2_lt_to_annot = pd.concat([s3_lt_to_annot, s4_lt_to_annot], axis=1)

    # Group #2 | Merge Prokka & PGAP annotations
    g2_prokka_pgap = pd.concat([g2_filt_gpa_reordered, g2_lt_to_annot], axis=1)

    # -------------------------------------------------------------------------------------------------
    # 3. Match other group Blastp results to validate Only-core-NF
    # The only-cores are compared once more to the annotations of other group strains through Blastp
    # -------------------------------------------------------------------------------------------------
    # Group #1 | G1 Only-core to G2 strains Blastp results
    ## G2-S3
    g1_oc_to_s3 = pd.read_csv(workspace + g1s3, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits
    g1_oc_to_s3_best = best_hit_oc_nf(df=g1_oc_to_s3, gpa_list=s1_gpa_list, 
                                      pident_cut=pident_cut_nf, eval_cut=evalue_cut_nf)

    ## Get only the columns needed & Rename columns
    header_change = ["sseqid", "pident", "evalue", "qcovs", "length"]
    g1_oc_to_s3_best = g1_oc_to_s3_best[header_change]
    rename = ["ID " + "[" + group_info2[0] + "]"]
    rename = rename + ["pident", "evalue", "qcovs", "length"]
    g1_oc_to_s3_best = g1_oc_to_s3_best.set_axis(rename, axis=1)

    ## G2-S4
    g1_oc_to_s4 = pd.read_csv(workspace + g1s4, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits
    g1_oc_to_s4_best = best_hit_oc_nf(df=g1_oc_to_s4, gpa_list=s1_gpa_list,
                                      pident_cut=pident_cut_nf, eval_cut=evalue_cut_nf)

    ## Get only the columns needed & Rename columns
    g1_oc_to_s4_best = g1_oc_to_s4_best[header_change]
    rename = ["ID " + "[" + group_info2[1] + "]"]
    rename = rename + ["pident", "evalue", "qcovs", "length"]
    g1_oc_to_s4_best = g1_oc_to_s4_best.set_axis(rename, axis=1)

    # Group #1 | Merge OC hit results
    g1_oc_to_other = pd.concat([g1_oc_to_s3_best, g1_oc_to_s4_best], axis=1)

    # 

    # Group #2 | G2 Only-core to G1 strains Blastp results
    ## G1-S1
    g2_oc_to_s1 = pd.read_csv(workspace + g2s1, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits
    g2_oc_to_s1_best = best_hit_oc_nf(df=g2_oc_to_s1, gpa_list=s3_gpa_list,
                                      pident_cut=pident_cut_nf, eval_cut=evalue_cut_nf)

    ## Get only the columns needed & Rename columns
    header_change = ["sseqid", "pident", "evalue", "qcovs", "length"]
    g2_oc_to_s1_best = g2_oc_to_s1_best[header_change]
    rename = ["ID " + "[" + group_info1[0] + "]"]
    rename = rename + ["pident", "evalue", "qcovs", "length"]
    g2_oc_to_s1_best = g2_oc_to_s1_best.set_axis(rename, axis=1)

    ## G1-S2
    g2_oc_to_s2 = pd.read_csv(workspace + g2s2, sep="\t", header=None, names=blastp_header)

    ## Parse only the best hits
    g2_oc_to_s2_best = best_hit_oc_nf(df=g2_oc_to_s2, gpa_list=s3_gpa_list,
                                      pident_cut=pident_cut_nf, eval_cut=evalue_cut_nf)

    ## Get only the columns needed & Rename columns
    g2_oc_to_s2_best = g2_oc_to_s2_best[header_change]
    rename = ["ID " + "[" + group_info1[1] + "]"]
    rename = rename + ["pident", "evalue", "qcovs", "length"]
    g2_oc_to_s2_best = g2_oc_to_s2_best.set_axis(rename, axis=1)

    # Group #2 | Merge OC hit results
    g2_oc_to_other = pd.concat([g2_oc_to_s1_best, g2_oc_to_s2_best], axis=1)

    # -------------------------------------------------------------------------------------------------
    # 4. Merge PGAP matching data and OC Blastp data
    # Final step to merge the Prokka-PGAP matching dataframe and OC-to-other group Blastp results
    # -------------------------------------------------------------------------------------------------
    # Group #1
    g1_final = pd.concat([g1_prokka_pgap, g1_oc_to_other], axis=1)

    # Group #2 
    g2_final = pd.concat([g2_prokka_pgap, g2_oc_to_other], axis=1)

    # -------------------------------------------------------------------------------------------------
    # 5. Save the final data
    # -------------------------------------------------------------------------------------------------
    # Group #1
    g1_final.to_csv(workspace + save_path + save_name1, sep=",", index=False)

    # Group #2 
    g2_final.to_csv(workspace + save_path + save_name2, sep=",", index=False)

    # -------------------------------------------------------------------------------------------------
    # 6. Print final message :)
    # -------------------------------------------------------------------------------------------------
    print("All tasks completed without error !")


if __name__ == "__main__":
    main()