# Optimisations
## Bowtie2
- Using gzipped files from the previous steps
- splitting the input and running in parallel before merging again. e.g. 
```
split -l 4000000 reads.fq split_reads_
parallel bowtie2 -x index -U {} -S {.}.sam ::: split_reads_*
```
- Setting more threads, somewhere between 10-16. Possibly allocate 2 more CPUs then threads according to biowulf documentation

## Samtools sort
- Bowtie2 pipes into `samtools sort`, but the sort doesn't use any multithreading. Can easily implement that with `-@ 8`
- Allocate more memory to the sorting process using the -m option. This option specifies the maximum memory per thread. e.g. `-m 4G`
- Might be able to write to /tmp for faster i/o e.g. `-T /tmp/$SLURM_JOB_ID`

## Recording stats
- this is seemingly done three times, once before replacing the read groups, and once after. Which of these steps are needed? presumably you just need that run once?

## GATK Haplotypecaller
- This natively supports a scatter-gather approach with the -L flag, which accepts a list of intervals. These need to be created in a previous step. The outputs need to be merged at the end as well in another seperate step. There is GATK documentation for all of this, although it can be hard to digest.

## GATK MarkDuplicates
- this is less of a priority
- tends to be i/o bound, so fast disk is ideal for the /tmp dir. We might be able to set this to the system /tmp if it's not too big. 
- Online I saw that setting `SORTING_COLLECTION_SIZE_RATIO` to something like 0.15 or 0.1 might help speed things up
- Could also try `MarkDuplicatesSpark` instead, which is meant to offer a 15% speedup
	- The tool is optimized to run on queryname-grouped alignments. If provided coordinate-sorted alignments, the tool will spend additional time first queryname sorting the reads internally. This can result in the tool being up to 2x slower processing under some circumstances.
	- Will scale linearly to upwards of 16 cores
	- MarkDuplicates run locally specifying the core input. Note if 'spark.executor.cores' is unset, Spark will use all available cores on the machine. Example: 
		```
		gatk MarkDuplicatesSpark \
            -I input.bam \
            -O marked_duplicates.bam \
            -M marked_dup_metrics.txt \
            --conf 'spark.executor.cores=5'
		```
-can avoid using samtools to create another index with ` -CREATE_INDEX true`
- Could also setup java optimisations
## AddOrReplaceReadGroups
- Is this step needed? Doesn't seem correct to overwrite all the fields with the alpaca ID
- Can create index at the same time with `CREATE_INDEX=true`and therefore skip the second samtools index process
- this job doesn't specifically need to have an indexed bam as input, could skip in indexing step if not too slow
- What about samtools flagstat and samtools idxstats- do we need that? 
