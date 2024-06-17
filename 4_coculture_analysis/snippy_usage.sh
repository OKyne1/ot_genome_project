# Code used in snippy work

singularity shell --bind /well/moru-batty/projects/ot/assembly/240530_ot_samples/coculturing:/well/moru-batty/projects/ot/assembly/240530_ot_samples/coculturing snippy.sif
## The --bind part is used to mount the required files, even if they are in the same location as the container


snippy --outdir 'snippy_multi_18-gil' --se 'gilliam/ot_and_unmapped_gilliam.fastq' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_19-gil' --se 'kato/ot_and_unmapped_kato.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_20-gil' --se 'a1_20/ot_and_unmapped_a1_20.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_21-gil' --se 'a2_21/ot_and_unmapped_a2_21.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_22-gil' --se 'c1_22/ot_and_unmapped_c1_22.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_23-gil' --se 'c2_23/ot_and_unmapped_c2_23.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05
snippy --outdir 'snippy_multi_karp-gil' --se '240409_combined_karp.fastq.gz' --ref data/Ot_gilliam_SAMEA104570318.fasta --prefix core --minfrac 0.05

snippy-core --ref 'snippy_multi_18-gil/ref.fa' snippy_multi_18-gil snippy_multi_20-gil snippy_multi_21-gil snippy_multi_22-gil snippy_multi_23-gil snippy_multi_19-gil snippy_multi_karp-gil
