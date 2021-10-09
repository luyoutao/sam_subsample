## Summary:
Random sample --num reads (SE) or read pairs (PE) from BAM or SAM. It's supposed to overcome the current limitation of `samtools view -s seed.perc` which cannot
output a fixed number of reads (e.g. https://www.biostars.org/p/145820/, https://bioinformatics.stackexchange.com/questions/402/how-can-i-downsample-a-bam-file-while-keeping-both-reads-in-pairs/406).

## Usage:
    sam_subsample --inFile input.[bam|sam] --outFile output.bam [--num 5000] [--seed 43] [--help] [--version] [--debug info]

## Options:
    -i, --infile FILE   input BAM/SAM, must be name sorted (@HD SO:queryname)
    -o, --outfile FILE  output BAM
    -n, --num INTEGER   number of reads (read pairs if PE) to downsample
                        (default: 5000)
    -s, --seed INTEGER  seed (default: None)
        --level         level of debugging info, choose from 'error', 'warn',
                        'info', 'debug', 'trace'
    -h, --help          print usage
    -v, --version       print version
    
