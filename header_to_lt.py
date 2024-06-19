# Date: 041724
import argparse
import os


def header_change(input, output):
    count = 0
    with open(input, "r") as infile, open(output, "w") as outfile:
    # Iterate through each line in the input file
        for line in infile:
            # Check if the line is a header line (starts with ">")
            if line.startswith(">"):
                c = line.count('>') - 1
                count += c
                if c > 0:
                    print(line)
                # Extract the locus tag from the header line
                locus_tag = line.split("[locus_tag=")[1].split("]")[0]
                # Write the new header line with only the locus tag
                outfile.write(">" + locus_tag + "\n")
            else:
                # Write sequence lines as they are
                outfile.write(line)
    return count


def main():    
    # Argument parser
    parser = argparse.ArgumentParser(description='Change the header of PGAP CDS protein FASTA to Locus Tags')
    parser.add_argument('--i', required=True, help='Path to raw FASTA file')
    parser.add_argument('--save_path', required=False, default='./', help='Path to save new FASTA file')
    parser.add_argument('--save_name', required=False, default='lt.fasta', help='Name of new FASTA file') 
    args = parser.parse_args()

    # Define argument variables
    cds = args.i
    save_path = args.save_path
    save_name = args.save_name

    output = os.path.join(save_path, save_name)
    c = header_change(cds, output)
    if c == 0:
        print("Successfully Converted!")

if __name__ == "__main__":
    main()
