import sys
import pysam
# from collections import Counter

def return_uniquely_STAR_mapped_reads(bam_in):
    out_count = 0
    for read in bam_in:
        if read.mapping_quality == 255:
            bam_out.write(read)
            out_count += 1
    print('reads_written =', out_count)
    bam_out.close()
    bam_in.close()

def return_under4loc_multi_STAR_mapped_reads(bam_in):
    out_count = 0
    for read in bam_in:
        if read.tags[0][1] in list(range(2,5)):
            bam_out.write(read)
            out_count += 1
    print('reads_written =', out_count)
    bam_out.close()
    bam_in.close()

def return_under10loc_multi_STAR_mapped_reads(bam_in):
    out_count = 0
    for read in bam_in:
        if read.tags[0][1] in list(range(2,11)):
            bam_out.write(read)
            out_count += 1
    print('reads_written =', out_count)
    bam_out.close()
    bam_in.close()

def return_ALL_multi_STAR_mapped_reads(bam_in):
    """
     Note that NH:i: tag in STAR will still report the actual number of loci that
    the reads map to. NH is tags[0].
    """

    out_count = 0
    for read in bam_in:
        # if read.tags[0][1] in list(range(2,21)):
        if read.tags[0][1] > 1:
            bam_out.write(read)
            out_count += 1
    print('reads_written =', out_count)
    bam_out.close()
    bam_in.close()

def return_ALL_multi_Bowtie_mapped_reads(bam_in):
    out_count = 0
    for read in bam_in:
        if read.tags[3][1] != 255:
            bam_out.write(read)
            out_count += 1
    print('reads_written =', out_count)
    bam_out.close()
    bam_in.close()


def return_mm_info_STAR(bam_in):
    counts = []
    for read in bam_in:
        NH = read.tags[0][1] # tags[0][1] is NH value.
        if NH > 1:
            counts.append(NH)

    counter_list = Counter(counts).most_common()

    print("Amount of reads that map multiple times: ", sum([item[1]/item[0] for item in counter_list]))

    print("Distribution of alignement count: ", counter_list)

    print("Total amount of alignements of ambiguously mapped reads: ", sum([item[1] for item in counter_list]))


if __name__ == '__main__':

    bam_in = pysam.Samfile(sys.argv[1], 'rb')
    bam_out = pysam.Samfile(sys.argv[2], 'wb', template=bam_in)

    # Amount of multi-mapped reads to output. 'under10' if used default --outFilterMultimapNmax with STAR
    # 'all' if used --outFilterMultimapNmax 20 in STAR
    amount = sys.argv[3]

    if amount == 'under4':
        return_under4loc_multi_STAR_mapped_reads(bam_in)
    if amount == 'under10':
        return_under10loc_multi_STAR_mapped_reads(bam_in)
    elif amount == 'all':
        return_ALL_multi_STAR_mapped_reads(bam_in)
