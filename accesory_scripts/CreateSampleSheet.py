#!/usr/bin/env python

"""
Create a samplesheet for the mgap pipeline from a directory of sequencing files.

Supports:
- Illumina paired-end reads (R1/R2)
- Oxford Nanopore long reads (single fastq)
- Pre-assembled contigs (fasta files)
"""

import argparse
import csv
import logging
import re
import sys
from pathlib import Path
from collections import defaultdict

logger = logging.getLogger()


class SampleSheetGenerator:
    """Generate samplesheet from a directory of sequencing files."""

    ILLUMINA_PATTERNS = [
        r"(.+?)[-_]S\d+[-_]",  # sample_S###_ (Illumina standard)
        r"(.+?)[-_]R[12][-_.]",  # sample_R1_ or sample-R1.
        r"(.+?)[-_]R[12]$",      # sample_R1 or sample-R1
        r"(.+?)[-_][12][-_.]",   # sample_1_ or sample-1.
        r"(.+?)[-_][12]$",       # sample_1 or sample-1
    ]

    ILLUMINA_R1_PATTERNS = [
        r"R1[-_.]",
        r"R1\.f(ast)?q",
        r"_1[-_.]",
        r"_1\.f(ast)?q",
    ]

    ILLUMINA_R2_PATTERNS = [
        r"R2[-_.]",
        r"R2\.f(ast)?q",
        r"_2[-_.]",
        r"_2\.f(ast)?q",
    ]

    FASTQ_EXTENSIONS = [".fastq.gz", ".fq.gz", ".fastq", ".fq"]
    FASTA_EXTENSIONS = [".fasta", ".fa", ".fna", ".fasta.gz", ".fa.gz", ".fna.gz"]

    def __init__(self, input_dir, output_file, data_type="auto"):
        """
        Initialize the samplesheet generator.

        Args:
            input_dir (Path): Directory containing sequencing files
            output_file (Path): Output samplesheet CSV file
            data_type (str): Type of data - 'illumina', 'ont', 'contig', or 'auto'
        """
        self.input_dir = Path(input_dir)
        self.output_file = Path(output_file)
        self.data_type = data_type
        self.samples = []

    def extract_sample_name(self, filename, is_illumina=False):
        """
        Extract sample name from filename.

        Args:
            filename (str): The filename to parse
            is_illumina (bool): Whether this is Illumina paired-end data

        Returns:
            str: The extracted sample name
        """
        # Remove extensions
        name = filename
        for ext in self.FASTQ_EXTENSIONS + self.FASTA_EXTENSIONS:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break

        if is_illumina:
            # Try to extract sample name before R1/R2 or _1/_2
            for pattern in self.ILLUMINA_PATTERNS:
                match = re.search(pattern, name)
                if match:
                    return match.group(1)

        # If no pattern matched or not Illumina, use the whole name
        # Clean up common suffixes
        name = re.sub(r'[-_](R[12]|[12])$', '', name)
        
        # Replace spaces and special characters with underscores
        name = re.sub(r'[^\w\-.]', '_', name)
        
        return name

    def is_r1_file(self, filename):
        """Check if filename is an R1 file."""
        for pattern in self.ILLUMINA_R1_PATTERNS:
            if re.search(pattern, filename):
                return True
        return False

    def is_r2_file(self, filename):
        """Check if filename is an R2 file."""
        for pattern in self.ILLUMINA_R2_PATTERNS:
            if re.search(pattern, filename):
                return True
        return False

    def get_file_type(self, filepath):
        """
        Determine the type of sequencing file.

        Returns:
            str: 'fastq' or 'fasta' or None
        """
        filename = filepath.name
        
        for ext in self.FASTQ_EXTENSIONS:
            if filename.endswith(ext):
                return 'fastq'
        
        for ext in self.FASTA_EXTENSIONS:
            if filename.endswith(ext):
                return 'fasta'
        
        return None

    def detect_data_type(self, files):
        """
        Auto-detect the type of sequencing data.

        Args:
            files (list): List of file paths

        Returns:
            str: 'illumina', 'ont', or 'contig'
        """
        fastq_files = []
        fasta_files = []
        r1_files = []
        r2_files = []

        for f in files:
            file_type = self.get_file_type(f)
            if file_type == 'fastq':
                fastq_files.append(f)
                if self.is_r1_file(f.name):
                    r1_files.append(f)
                elif self.is_r2_file(f.name):
                    r2_files.append(f)
            elif file_type == 'fasta':
                fasta_files.append(f)

        # If we have fasta files, assume contigs
        if fasta_files:
            return 'contig'
        
        # If we have paired R1/R2 files, assume Illumina
        if r1_files and r2_files:
            return 'illumina'
        
        # If we only have fastq files without clear R1/R2 pairing, assume ONT
        if fastq_files:
            return 'ont'
        
        logger.warning("Could not auto-detect data type. Defaulting to 'ont'")
        return 'ont'

    def process_illumina(self, files):
        """Process Illumina paired-end reads."""
        # Group files by sample name
        samples_dict = defaultdict(lambda: {'R1': [], 'R2': []})
        
        for filepath in files:
            if self.get_file_type(filepath) != 'fastq':
                continue
            
            filename = filepath.name
            sample_name = self.extract_sample_name(filename, is_illumina=True)
            
            if self.is_r1_file(filename):
                samples_dict[sample_name]['R1'].append(str(filepath.absolute()))
            elif self.is_r2_file(filename):
                samples_dict[sample_name]['R2'].append(str(filepath.absolute()))
            else:
                logger.warning(f"Could not determine if {filename} is R1 or R2. Skipping.")
        
        # Create sample entries
        for sample_name, reads in sorted(samples_dict.items()):
            r1_files = sorted(reads['R1'])
            r2_files = sorted(reads['R2'])
            
            if not r1_files:
                logger.warning(f"No R1 files found for sample {sample_name}. Skipping.")
                continue
            
            # Match R1 and R2 files
            for i, r1 in enumerate(r1_files):
                r2 = r2_files[i] if i < len(r2_files) else ""
                
                if not r2:
                    logger.warning(f"No matching R2 for {r1}. Adding as single-end.")
                
                self.samples.append({
                    'sample': sample_name,
                    'fastq_1': r1,
                    'fastq_2': r2
                })

    def process_ont(self, files):
        """Process Oxford Nanopore long reads."""
        for filepath in sorted(files):
            if self.get_file_type(filepath) != 'fastq':
                continue
            
            sample_name = self.extract_sample_name(filepath.name, is_illumina=False)
            
            self.samples.append({
                'sample': sample_name,
                'fastq_1': str(filepath.absolute()),
                'fastq_2': ''
            })

    def process_contigs(self, files):
        """Process pre-assembled contigs."""
        for filepath in sorted(files):
            if self.get_file_type(filepath) != 'fasta':
                continue
            
            sample_name = self.extract_sample_name(filepath.name, is_illumina=False)
            
            self.samples.append({
                'sample': sample_name,
                'fasta': str(filepath.absolute())
            })

    def generate(self):
        """Generate the samplesheet."""
        if not self.input_dir.is_dir():
            logger.error(f"Input directory {self.input_dir} does not exist!")
            sys.exit(1)
        
        # Get all files in the directory
        all_files = [f for f in self.input_dir.iterdir() if f.is_file()]
        
        if not all_files:
            logger.error(f"No files found in {self.input_dir}")
            sys.exit(1)
        
        # Auto-detect data type if needed
        if self.data_type == 'auto':
            self.data_type = self.detect_data_type(all_files)
            logger.info(f"Auto-detected data type: {self.data_type}")
        
        # Process files based on data type
        if self.data_type == 'illumina':
            self.process_illumina(all_files)
        elif self.data_type == 'ont':
            self.process_ont(all_files)
        elif self.data_type == 'contig':
            self.process_contigs(all_files)
        else:
            logger.error(f"Unknown data type: {self.data_type}")
            sys.exit(1)
        
        if not self.samples:
            logger.error("No valid samples found!")
            sys.exit(1)
        
        # Write samplesheet
        self.write_samplesheet()
        
        logger.info(f"Successfully created samplesheet with {len(self.samples)} entries")
        logger.info(f"Output written to: {self.output_file}")

    def write_samplesheet(self):
        """Write the samplesheet to CSV file."""
        self.output_file.parent.mkdir(parents=True, exist_ok=True)
        
        # Use different columns for contig mode
        if self.data_type == 'contig':
            fieldnames = ['sample', 'fasta']
        else:
            fieldnames = ['sample', 'fastq_1', 'fastq_2']
        
        with self.output_file.open('w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(self.samples)


def parse_args(argv=None):
    """Define and parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate a samplesheet for the mgap pipeline from a directory of sequencing files.",
        epilog="Example: python CreateSampleSheet.py /path/to/fastq_dir samplesheet.csv --type illumina"
    )
    parser.add_argument(
        "input_dir",
        type=Path,
        help="Directory containing sequencing files (FASTQ or FASTA)"
    )
    parser.add_argument(
        "output_file",
        type=Path,
        help="Output samplesheet CSV file"
    )
    parser.add_argument(
        "-t", "--type",
        dest="data_type",
        choices=['auto', 'illumina', 'ont', 'contig'],
        default='auto',
        help="Type of sequencing data (default: auto-detect)"
    )
    parser.add_argument(
        "-l", "--log-level",
        choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG'],
        default='INFO',
        help="Logging level (default: INFO)"
    )
    
    return parser.parse_args(argv)


def main(argv=None):
    """Main entry point."""
    args = parse_args(argv)
    
    logging.basicConfig(
        level=args.log_level,
        format='[%(levelname)s] %(message)s'
    )
    
    generator = SampleSheetGenerator(
        input_dir=args.input_dir,
        output_file=args.output_file,
        data_type=args.data_type
    )
    
    generator.generate()


if __name__ == "__main__":
    sys.exit(main())
