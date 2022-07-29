### Background

tRNA quantification is hard. cDNA synthesis is hampered by:

  * pervasive RT blocks by modified (mostly methylated) nucleosides ... this leads to premature RT stops at modified sites... resulting in very short reads.
  * extensive similarity between tRNA genes ... eg. tRNA genes coding for 2 different codons can have only a single nucleotide difference. ... this make the _alignment_ step of analysis pipeline extremely difficult.
  
The composition of cellular tRNA pools is critical for efficient mRNA decoding and proteome integrity. Defective tRNA biogenesis is linked to neurological disorders and cancer.

<u>AlkB</u>: AlkB is a protein found in E. coli which removes alkylation groups from ssDNA.

<u>Isodecoder</u>:  tRNA genes that share the same anticodon but have different body sequences have been termed as “isodecoders”, implying that they all read the same codon in translation.

<u>tRNA human genes</u>: Facts:
  
  * Number of different tRNA anticodon families: 57
  * Number of unique tRNA molecules produced (all cell types included): ~400

<u>K562</u>: Human CML (chronic myelogenous leukemia) cell line

### mimseq and its issues

[mimseq](https://www.sciencedirect.com/science/article/pii/S1097276521000484) is the analysis tool to be used for quantifying tRNA abundance for Wu lab samples (K562 cell line samples). Code can be found [here](https://github.com/nedialkova-lab/mim-tRNAseq)

#### "group argument must be None for now"

After installation, got this error:

```bash
% mimseq  --species Hsap  --cluster-id 0.95  --threads 15  --min-cov 0.0005  --max-mismatches 0.
1  --control-condition HEK293T  -n hg38_test  --out-dir hg38_HEK239vsK562  --max-multi 4 --rem
ap  --remap-mismatches 0.075  sampleData_HEKvsK562.txt
           _                 _   ____  _   _    _
 _ __ ___ (_)_ __ ___       | |_|  _ \| \ | |  / \   ___  ___  __ _
| '_ ` _ \| | '_ ` _ \ _____| __| |_) |  \| | / _ \ / __|/ _ \/ _` |
| | | | | | | | | | | |_____| |_|  _ <| |\  |/ ___ \\__ \  __/ (_| |
|_| |_| |_|_|_| |_| |_|      \__|_| \_\_| \_/_/   \_\___/\___|\__, |
                                                                 |_|

 Modification-induced misincorporation analysis of tRNA sequencing data

2022-06-27 15:51:26,542 [INFO ] mim-tRNAseq v1.1.6 run with command:
2022-06-27 15:51:26,542 [INFO ] /usr/local/bin/mimseq --species Hsap --cluster-id 0.95 --threads 15 --min-cov 0.0005 --max-mismatches 0.1 --control-condition HEK293T -n hg38_test --out-dir hg38_HEK239vsK562 --max-multi 4 --remap --remap-mismatches 0.075 sampleData_HEKvsK562.txt
2022-06-27 15:51:26,544 [INFO ]
...
...
+--------------------------------------------------------------------+
| Analysing misincorporations and stops to RT, and analysing 3' ends |
+--------------------------------------------------------------------+
2022-06-27 15:56:15,178 [INFO ] ** Discovering unannotated modifications for realignment **
Traceback (most recent call last):
  File "/usr/local/bin/mimseq", line 8, in <module>
    sys.exit(main())
  File "/usr/local/lib/python3.8/dist-packages/mimseq/mimseq.py", line 407, in main
    mimseq(args.trnas, args.trnaout, args.name, args.species, args.out, args.cluster, args.cluster_id, args.cov_diff, \
  File "/usr/local/lib/python3.8/dist-packages/mimseq/mimseq.py", line 132, in mimseq
    new_mods, new_Inosines, filtered_cov, filter_warning, unsplitCluster_lookup,readRef_unsplit_newNames = generateModsTable(coverageData, out, name, threads, min_cov, mismatch_dict, insert_dict, del_dict, cluster_dict, cca, remap, misinc_thresh, mod_lists, Inosine_lists, tRNA_dict, Inosine_clusters, unique_isodecoderMMs_new, splitBool_new, isodecoder_sizes, unsplitCluster_lookup, cluster)
  File "/usr/local/lib/python3.8/dist-packages/mimseq/mmQuant.py", line 634, in generateModsTable
    pool = MyPool(multi)
  File "/usr/lib/python3.8/multiprocessing/pool.py", line 212, in __init__
    self._repopulate_pool()
  File "/usr/lib/python3.8/multiprocessing/pool.py", line 303, in _repopulate_pool
    return self._repopulate_pool_static(self._ctx, self.Process,
  File "/usr/lib/python3.8/multiprocessing/pool.py", line 319, in _repopulate_pool_static
    w = Process(ctx, target=worker,
  File "/usr/lib/python3.8/multiprocessing/process.py", line 82, in __init__
    assert group is None, 'group argument must be None for now'
AssertionError: group argument must be None for now
```

This was resolved with specific version of python/conda/mamba. See [here](https://github.com/nedialkova-lab/mim-tRNAseq/issues/33) for details.

#### Keyerror

With the above issue fixed, another error was found

```bash
+-----------------------------------+
| Final deconvolution and filtering |
+-----------------------------------+
2022-07-10 12:04:32,132 [INFO ] 234 clusters filtered out according to minimum coverage threshold: 0.05% of total tRNA coverage.
2022-07-10 12:05:19,876 [INFO ] Output final tables, counts and new mods...
2022-07-10 12:05:23,815 [INFO ] ** Read counts per anticodon saved to Tet2KOvsWT/counts/Anticodon_counts_raw.txt **
2022-07-10 12:05:23,815 [INFO ] ** Read counts per isodecoder saved to Tet2KOvsWT/counts/Isodecoder_counts_raw.txt **
2022-07-10 12:05:23,815 [INFO ] 11 sequences could not be deconvoluted as >10% parent assigned-reads do not match parent sequence!
2022-07-10 12:05:23,815 [INFO ] Reasons for this include misaligment of reads to an incorrect cluster, or inaccurate aligment by GSNAP (e.g. at indels in reads) prohibiting correct devonvolution.
2022-07-10 12:05:23,815 [INFO ] 79 total unique sequences not deconvoluted due to mismatches at modified sites, insufficient coverage or read mismatches to parent
Traceback (most recent call last):
  File "/opt2/conda/envs/mimseq/bin/mimseq", line 10, in <module>
    sys.exit(main())
  File "/opt2/conda/envs/mimseq/lib/python3.7/site-packages/mimseq/mimseq.py", line 410, in main
    args.misinc_thresh, args.mito, args.pretrnas, args.local_mod, args.p_adj, args.sampledata)
  File "/opt2/conda/envs/mimseq/lib/python3.7/site-packages/mimseq/mimseq.py", line 161, in mimseq
    isodecoder_sizes = writeIsodecoderInfo(out, name, isodecoder_sizes,readRef_unsplit_newNames, tRNA_dict)
  File "/opt2/conda/envs/mimseq/lib/python3.7/site-packages/mimseq/splitClusters.py", line 477, in writeIsodecoderInfo
    del isodecoder_sizes[gene]
KeyError: 'Homo_sapiens_tRNA-Glu-TTC-1-1'
```

This issue was reported by other users as well and the author recommended using the latest updated version of `splitClusters.py`. See [here](https://github.com/nedialkova-lab/mim-tRNAseq/issues/34) for details.

