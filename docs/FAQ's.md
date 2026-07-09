<details>
<summary><b>Whats the difference between --read-wide and --loci-wide mode?</b></summary>

By default ATaRVa uses `--read-wise` mode which processes each read in BAM and stores required information, but in `--loci-wise` mode it processes TR loci based on the order given in the refernce BED file.
You should use the `--loci-wise` when your reference BED TR loci are sparsely scattered across genome, which can reduce the overall genotyping time. 

</details>