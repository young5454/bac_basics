# Date: 061324
import argparse
import os
import pandas as pd


def make_meta(input, output):
    meta, lt, ds = {}, [], []
    with open(input, "r") as cds_fasta, open(output, "w") as meta_csv:
        for line in cds_fasta:
            # Check if the line is a header line (starts with ">")
            if line.startswith(">"):
                # Extract the locus tag from the header line
                locus_tag = line.split("[locus_tag=")[1].split("]")[0]
                description = line.split("[protein=")[1].split("]")[0]
                lt.append(locus_tag)
                ds.append(description)
        # Save data into dictionary
        meta['LocusTag'] = lt
        meta['Description'] = ds
        # Change dictionary into dataframe
        meta_df = pd.DataFrame(meta)

    # Save metadata
    meta_df.to_csv(output, index = False)

    return meta_df


def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Make matched metadata of LocusTag : Description for NCBI PGAP annotation')
    parser.add_argument('--i', required=True, help='Path to NCBI PGAP-annotated CDS protein FASTA file')
    parser.add_argument('--save_path', required=False, default='./', help='Path to save new FASTA file')
    parser.add_argument('--save_name', required=False, default='pgap_meta.csv', help='Name of metadata file' ) 
    args = parser.parse_args()

    # Define argument variables
    cds = args.i
    save_path = args.save_path
    save_name = args.save_name

    output = os.path.join(save_path, save_name)
    meta_df = make_meta(cds, output)

if __name__ == "__main__":
    main()