#!/usr/bin/env nextflow

// =============================================================================
// bionets.nf  -  Step 4 of gwas-bionets: biological network (module) methods
// =============================================================================
// Runs network-based methods that combine gene-level scores (Step 3) with a
// reference interaction network to extract disease-associated modules. Each method
// is a separate Nextflow module in the same folder (called as ${baseDir}/<method>.nf)
// and is run once per split (or chunk) I = 1..K.
//
//   Method                 Module        Called from
//   Heinz                  heinz.nf      heinz_call
//   dmGWAS                 dmgwas.nf     dmgwas_call
//   SigMod                 sigmod.nf     sigmod_call
//   HotNet2                hotnet2.nf    hotnet2_call
//   LEAN                   lean.nf       lean_call
//   Hierarchical HotNet    hhotnet.nf    hhotnet_call (currently DISABLED in the workflow)
//
// Where the scores are read and the results written (<scores> = base name of --scores):
//   k = 1                : <out>/<scores>.tsv                          -> <out>/<method>
//   k > 1, --d_samp 1    : <out>/stable_consensus/data_splits/data_split_<I>/<scores>_split_<I>.tsv
//                          -> <out>/stable_consensus/data_splits/data_split_<I>/<method>
//   k > 1, --d_samp 0    : <out>/scores_chunks/scores_chunk_<I>/<scores>_chunk_<I>.tsv
//                          -> <out>/scores_chunks/scores_chunk_<I>/<method>
//
// Parameters (command line, as --name value)
//   General
//     --network   (required) Reference network: tab-separated list of interacting gene pairs
//                 (undirected, no self-loops or duplicated edges, gene symbols matching the scores).
//     --scores    (required) Gene scores file produced by magma_calc.nf (gene_ids_reformat),
//                 '<prefix>.scores.genes.out_converted.tsv'; its base name is used to build the
//                 per-split file names above.
//     --out       Output folder (default '.').
//     --k         Number of splits (data splits or score chunks); it must match the previous steps (default 1).
//     --d_samp    1 (default) if the k samples come from the data (genotype splits, Steps 1-3);
//                 0 if the gene scores file was divided into k chunks and each chunk is used as scores file.
//   Heinz
//     --heinz_fdr           False discovery rate used to score nodes from the Beta-uniform mixture
//                           model; lower values give smaller subnetworks (default 0.2).
//   dmGWAS
//     --dmgwas_d            Order of the neighbours searched when expanding a module (default 2;
//                           1 = direct neighbours only). Order 2 is only available in dmGWAS_1.0/2.X.
//     --dmgwas_r            Minimum relative score increment to add a neighbour to a module
//                           (default 0.1; larger values are stricter).
//   SigMod
//     --sigmod_path         Folder with SigMod's code (SigMod_v2, R functions).
//     --sigmod_lambdamax    Upper bound of the search interval for lambda (default 1)  } defined here
//     --sigmod_nmax         Maximum number of genes considered (default 300)           } but NOT forwarded
//     --sigmod_maxjump      Minimum size jump of the selection path (default 10)       } to sigmod.nf
//   HotNet2
//     --lfdr                Local FDR cutoff used to sparsify the scores: genes with a P-value >= cutoff get no heat.
//     --hotnet2_path        Folder with HotNet2's code.
//   LEAN
//     --reps                Number of random P-values drawn for each background distribution
//                           (n_reps in LEANR::run.lean, default 10000 there; 100000 here).
//   Hierarchical HotNet (inactive while hhotnet_call stays commented out in the workflow)
//     --hhotnet_path        Folder with the Hierarchical HotNet code.
//     --conn                Connectivity setting, forwarded to hhotnet.nf as --connectivity.
//     --perms               Number of permutations used to assess significance.
//     --beta                Restart probability of the random walk used for the similarity matrix.
//   --help                  Print the usage and exit.
//
// NOTE: every call below hard-codes '-profile cbio_cluster' for the child module, so the
//       profile given to bionets.nf is not propagated.

// General parameters
params.out = '.'             // --out: output folder
params.k = 1                 // --k: number of splits or chunks
params.d_samp = 1            // --d_samp: 1 = splits from the data, 0 = chunks of the scores file
// (no default: --network, --scores)
//params.profile = 'local'

// SigMod parameters (defined, but not forwarded to sigmod.nf: SigMod runs with its own defaults)
// (no default: --sigmod_path)
params.sigmod_lambdamax = 1   // upper bound of the lambda search interval
params.sigmod_nmax = 300      // maximum number of genes considered
params.sigmod_maxjump = 10    // minimum size jump in the selection path

// Heinz parameters
params.heinz_fdr = 0.2        // FDR for node scoring (lower = smaller subnetworks)

// dmGWAS parameters
params.dmgwas_r = 0.1         // minimum relative score increment to expand a module
params.dmgwas_d = 2           // order of the neighbours searched (1 or 2)

// HotNet2, LEAN and Hierarchical HotNet parameters have no default value here:
//   --hotnet2_path          (HotNet2)
//   --hhotnet_path  (Hierarchical HotNet, currently disabled)
params.lfdr = 0.125

// LEAN parameters
params.reps = 100000

// Hierarchical Hotnet parameters
params.conn = 1
params.perms = 100
params.beta = 0.5

// Help info ##########
params.help = null            // --help: print the usage and exit

/////////////////////////////////////
//  HEINZ METHOD
/////////////////////////////////////

// heinz_call: one task per split I. Runs heinz.nf with the split's scores, the network (--tab2),
// the FDR (--fdr), the split index (--i) and --d_samp; results in <...>/heinz.
process heinz_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val fdr
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/heinz.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --out ${out}${dir_splits}/heinz \
        --fdr ${fdr} \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  DMGWAS METHOD
/////////////////////////////////////

// dmgwas_call: one task per split I. Runs dmgwas.nf with the split's scores, the network (--tab2),
// d (--d) and r (--r); results in <...>/dmgwas.
process dmgwas_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val d
      val r
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/dmgwas.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --out ${out}${dir_splits}/dmgwas \
        --d ${d} \
        --r ${r} \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}


/////////////////////////////////////
//  SIGMOD METHOD
/////////////////////////////////////

// sigmod_call: one task per split I. Runs sigmod.nf with the split's scores, the network (--network)
// and the SigMod code folder (--sigmod_path); lambda_max, nmax and maxjump are not passed on.
// Results in <...>/sigmod.
process sigmod_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val sigmod_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/sigmod.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --sigmod_path ${sigmod_path} \
        --out ${out}${dir_splits}/sigmod \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  HOTNET2 METHOD
/////////////////////////////////////

// hotnet2_call: one task per split I. Runs hotnet2.nf with the split's scores, the network (--tab2),
// the local FDR cutoff (--lfdr_cutoff) and the HotNet2 code folder; results in <...>/hotnet2.
process hotnet2_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val lfdr
      val hotnet2_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/hotnet2.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --tab2 ${network} \
        --lfdr_cutoff ${lfdr} \
        --hotnet2_path ${hotnet2_path} \
        --out ${out}${dir_splits}/hotnet2 \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  HIERARCHICAL HOTNET METHOD
/////////////////////////////////////

// hhotnet_call: one task per split I. Runs hhotnet.nf with the split's scores, the network, connectivity,
// permutations, beta and the Hierarchical HotNet code folder; results in <...>/hhotnet.
// Not executed at the moment (see the workflow).
process hhotnet_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val conn
      val perms
      val beta
      val hhotnet_path
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/hhotnet.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --connectivity ${conn} \
        --permutations ${perms} \
        --beta ${beta} \
        --hhotnet_path ${hhotnet_path} \
        --out ${out}${dir_splits}/hhotnet \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}

/////////////////////////////////////
//  LEAN METHOD
/////////////////////////////////////

// lean_call: one task per split I. Runs lean.nf with the split's scores, the network and the number of
// repetitions (--reps) of the LEAN background distributions; results in <...>/lean.
process lean_call {
    input:
      file scores
      val network
      val k
      val d_samp
      path out
      val reps
      each I

    script:
      def split_suffix = (k > 1 && d_samp != 0) ? "_split_${I}" : (k > 1 && d_samp == 0) ? "_chunk_${I}" : ""
      def dir_splits = (k > 1 && d_samp != 0) ? "/stable_consensus/data_splits/data_split_${I}" : (k > 1 && d_samp == 0) ? "/scores_chunks/scores_chunk_${I}" : ""

      """
      ${baseDir}/lean.nf \
        --scores ${out}${dir_splits}/${scores.baseName}${split_suffix}.tsv \
        --network ${network} \
        --reps ${reps} \
        --out ${out}${dir_splits}/lean \
        --i ${I} \
        --d_samp ${d_samp} \
        -profile cbio_cluster
      """
}


// Workflow
workflow {
  log.info ""
  log.info "========================================================"
  log.info "|          [gwas-bionets] - bionets.nf          |"
  log.info "========================================================"
  log.info ""
  log.info "### Script to call different biological network methods ###"
  log.info ""
  log.info ""
  log.info ""
  log.info ""
  log.info "--------------------------------------------------------"
  log.info "This program comes with NO WARRANTY"
  log.info "It is free software, see LICENSE for details about"
  log.info "redistribution and contribution."
  log.info "--------------------------------------------------------"
  log.info ""

  if (params.help) {
    log.info "Usage : bionets.nf --network <network_file> --scores <scores_file> --k <knumber> --d_samp <data_sampling> --out <filename> \\"
    log.info "                   [--heinz_fdr <fdr>] [--dmgwas_d <d>] [--dmgwas_r <r>] [--lfdr <lfdr_cutoff>] [--reps <n>] \\"
    log.info "                   [--sigmod_path <path>] [--hotnet2_path <path>]"
    log.info ""
    log.info "  --network       (required) network of reference for the methods to operate and build a network according to gene"
    log.info "                  scores (p-values): tab-separated list of interacting gene pairs (undirected, without self-loops)."
    log.info "  --scores        (required) file for the formatted gene scores, produced by the 'gene_ids_reformat' process of the"
    log.info "                  magma_calc.nf file. Its base name is used to find the file of each split or chunk."
    log.info "  --out           path/to/folder where the results will be saved (default '.')."
    log.info "  --k             number of times data was split, it must be in accordance to the previous k for splitting the data."
    log.info "  --d_samp        1 or 0 (default 1) to specify where it was sampled from. For instance, you sample the data in different"
    log.info "                  splits (1) or one can divide the gene scores file in different chunks and use each chunk as your"
    log.info "                  gene scores (or p-values) file (0)."
    log.info "  --heinz_fdr     false discovery rate parameter to control the resultant subnetwork size in the Heinz method"
    log.info "                  (default 0.2; lower values give smaller subnetworks)."
    log.info "  --dmgwas_d      dmGWAS: order of the neighbour genes searched (default 2; dmGWAS_3.0 only uses 1)."
    log.info "  --dmgwas_r      dmGWAS: minimum relative score increment to add a neighbour to a module (default 0.1)."
    log.info "  --sigmod_path   SigMod: path to SigMod files with the method's internal code. Its parameters lambda_max, nmax and"
    log.info "                  maxjump (--sigmod_lambdamax, --sigmod_nmax, --sigmod_maxjump) are not forwarded: SigMod defaults apply."
    log.info "  --lfdr          local fdr cutoff of use in HotNet2 method to sparsify the scores (genes with P-value >= cutoff get no heat)."
    log.info "  --hotnet2_path  HotNet2: path to HotNet2 files with the method's internal code."
    log.info "  --reps          LEAN: number of random p-values drawn for each background distribution."
    log.info "  --hhotnet_path, --conn, --perms, --beta"
    log.info "                  Hierarchical HotNet: code folder, connectivity, number of permutations and restart probability."
    log.info "                  Currently inactive (the call is commented out in the workflow)."
    log.info "  --help          show this message and exit."
    log.info ""
    log.info "Example :"
    log.info "bionets.nf \\"
    log.info "  --network           path/to/my_network \\"
    log.info "  --k                 5 \\"
    log.info "  --d_samp            1 \\"
    log.info "  --scores            path/to/my_scores \\"
    log.info "  --heinz_fdr         0.5 \\"
    log.info "  --lfdr              0.125 \\"
    log.info "  --reps              10000 \\"
    log.info "  --out               path/to/my_output \\"
    log.info "  -profile            my_cluster \\"

    exit 0
  }

  // Define inputs
  def scores = file(params.scores)
  def network = file(params.network)
  def sigmod_path = params.sigmod_path
  def hotnet2_path = params.hotnet2_path
  def hhotnet_path = params.hhotnet_path
  def K = params.k as Integer
  def d_samp = params.d_samp as Integer
  def fdr = params.heinz_fdr as Double
  def lfdr = params.lfdr as Double
  def conn =  params.conn as Integer
  def perms = params.perms as Integer
  def beta = params.beta as Double
  def D = params.dmgwas_d as Integer
  def R = params.dmgwas_r as Double
  def lambdamax = params.sigmod_lambdamax as Integer
  def nmax = params.sigmod_nmax as Integer
  def maxjump = params.sigmod_maxjump as Integer
  def reps = params.reps as Integer

  // println "The value of scores is: ${params.scores}"
  // println "The value of out is: ${params.out}"
  // println "The value of scores baseName is: ${scores.baseName}"

  // Run processes inline
  heinz_call(scores, network, K, d_samp, params.out, fdr, 1..K)
  dmgwas_call(scores, network, K, d_samp, params.out, D, R, 1..K)
  sigmod_call(scores, network, K, d_samp, params.out, sigmod_path, 1..K)
  hotnet2_call(scores, network, K, d_samp, params.out, lfdr, hotnet2_path, 1..K)
  lean_call(scores, network, K, d_samp, params.out, reps, 1..K)
  // NOTE: Hierarchical HotNet is disabled (commented out); re-enable both
  // to run it. The Heinz, dmGWAS, SigMod, HotNet2 and LEAN calls above run for every split.
  // hhotnet_call(scores, network, K, d_samp, params.out, conn, perms, beta, hhotnet_path, 1..K)
  
}

