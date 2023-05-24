import os
import argparse

# Read arguments
parser = argparse.ArgumentParser(description="Script to generate a samplesheet for the pipeline")
parser.add_argument("-i", "--input_folder", help="Input folder with the fastq files", required=True)
parser.add_argument("-o", "--output_file", help="Output file with the samplesheet", required=True)

args = parser.parse_args()

# Get the fastq files
r1 = sorted([f for f in os.listdir(args.input_folder) if f.endswith("R1_001.fastq.gz")])
r2 = sorted([f for f in os.listdir(args.input_folder) if f.endswith("R2_001.fastq.gz")])


sample_name = [f.split("_")[0] for f in r1]

# Write the samplesheet
with open(args.output_file, "w") as f:
    f.write("sample,fastq_1,fastq_2\n")
    for i in range(len(r1)):
        f.write("{},{},{}\n".format(sample_name[i], os.path.join(args.input_folder, r1[i]), os.path.join(args.input_folder, r2[i])))


