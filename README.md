## Summary:
Random sample --num reads (SE) or read pairs (PE) from BAM or SAM

## Usage:
sam_subsample --inFile input.[bam|sam] --outFile output.bam [--num 5000] [--seed 43] [--help] [--version] [--debug info]

## Options:
    -i, --infile FILE   input BAM/SAM, queryname sorted
    -o, --outfile FILE  output BAM
    -n, --num INTEGER   number of reads (read pairs if PE) to downsample
                        (default: 5000)
    -s, --seed INTEGER  seed (default: None)
        --level         level of debugging info, choose from 'error', 'warn',
                        'info', 'debug', 'trace'
    -h, --help          print usage
    -v, --version       print version
    
