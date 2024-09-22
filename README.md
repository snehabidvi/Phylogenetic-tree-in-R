# Phylogenetic tree in R

A phylogenetic tree, also called phylogeny, is a diagram depicting the lines of evolutionary descent of different species, organisms, or genes from a common ancestor. Phylogenies are useful for organizing and visualizing knowledge of biological diversity, for structuring classifications, and for providing insight into events that occurred during evolution. As these trees show descent from a common ancestor and because much of the strongest evidence for evolution comes in the form of common ancestry, it is important to understand phylogenies for fully appreciating the overwhelming evidence that supports the theory of evolution.

This workflow is based on R programming, requiring packages such as Biostrings, DECIPHER, phangorn, ggtree and ggplot2

If these packages are not installed, then run the following commands to first install the packages.




## To install Bioconductor suite of packages, run -

```{r}
#if (!require("BiocManager", quietly = TRUE))
    #install.packages("BiocManager")
#BiocManager::install(version = "3.19")

```


## To install Biostrings, ggtree and DECIPHER packages which are the part of Bioconductor library, run -

```{r}
#BiocManager::install(c("Biostrings","ggtree" ,"DECIPHER"))

```


The input file of sequences must be in the fasta format. The sequence file should be saved in particular directory and the working directory should be specified correctly to import that file in R. The sequences should be named with the organisms' names so that those will be represented at the nodes of phylogenetic tree. 




## Reading XStringSet object from a file

```{r}
#readDNAStringSet(filepath, format="fasta",
               #nrec=-1L, skip=0L, seek.first.rec=FALSE,
               #use.names=TRUE, with.qualities=FALSE)

```

Arguments -

    1. filepath	:- A character vector (of arbitrary length when reading, of length 1 when writing) containing the path(s) to the file(s) to read or write. Reading files in gzip format (which usually have the '.gz' extension) is supported.

    Note that special values like "" or"|cmd" (typically supported by other I/O functions in R) are not supported here.

    Also filepath cannot be a standard connection. However filepath can be an object as returned by open_input_files. This object can be used to read files by chunks. 

    2. format :- Either "fasta" (the default) or "fastq".

    3. nrec :- Single integer. The maximum of number of records to read in. Negative values are ignored.

    4. skip :- Single non-negative integer. The number of records of the data file(s) to skip before beginning to read in records.

    5. seek.first.rec :- TRUE or FALSE (the default). If TRUE, then the reading function starts by setting the file position indicator at the beginning of the first line in the file that looks like the beginning of a FASTA (if format is "fasta") or FASTQ (if format is "fastq") record. More precisely this is the first line in the file that starts with a '>' (for FASTA) or a '@' (for FASTQ). An error is raised if no such line is found.

    Normal parsing then starts from there, and everything happens like if the file actually started there. In particular it will be an error if this first record is not a valid FASTA or FASTQ record.

    Using seek.first.rec=TRUE is useful for example to parse GFF3 files with embedded FASTA data.

    6. use.names	:- TRUE (the default) or FALSE. If TRUE, then the returned vector is named. For FASTA the names are taken from the record description lines. For FASTQ they are taken from the record sequence ids. Dropping the names with use.names=FALSE can help reduce memory footprint e.g. for a FASTQ file containing millions of reads.

    7. with.qualities :- TRUE or FALSE (the default). This argument is only supported when reading a FASTQ file. If TRUE, then the quality strings are also read and returned in the qualities metadata column of the returned DNAStringSet object. Note that by default the quality strings are ignored. This helps reduce memory footprint if the FASTQ file contains millions of reads.

    Note that using readQualityScaledDNAStringSet() is the preferred way to load a set of DNA sequences and their qualities from a FASTQ file into Bioconductor. Its main advantage is that it will return a QualityScaledDNAStringSet object instead of a DNAStringSet object, which makes handling of the qualities more convenient and less error prone. 

    8. seqtype	:- A single string specifying the type of sequences contained in the FASTA file(s). Supported sequence types:

    "B" for anything i.e. any letter is a valid one-letter sequence code.

    "DNA" for DNA sequences i.e. only letters in DNA_ALPHABET (case ignored) are valid one-letter sequence codes.

    "RNA" for RNA sequences i.e. only letters in RNA_ALPHABET (case ignored) are valid one-letter sequence codes.

    "AA" for Amino Acid sequences i.e. only letters in AA_ALPHABET (case ignored) are valid one-letter sequence codes.

    Invalid one-letter sequence codes are ignored with a warning.

    9. x	:- For writeXStringSet, the object to write to file.

    For saveXStringSet, the object to serialize.

    10. append :- TRUE or FALSE. If TRUE output will be appended to file; otherwise, it will overwrite the contents of file. See ?cat for the details.

    11. compress :- Like for the save function in base R, must be TRUE or FALSE (the default), or a single string specifying whether writing to the file is to use compression. The only type of compression supported at the moment is "gzip".

    Passing TRUE is equivalent to passing "gzip".

    12. compression_level	:- Not implemented yet.

    13. ...	:- Further format-specific arguments.

    If format="fasta", the width argument can be used to specify the maximum number of letters per line of sequence. width must be a single integer.

    If format="fastq", the qualities argument can be used to specify the quality strings. qualities must be a BStringSet object. If the argument is omitted, then the quality strings are taken from the qualities metadata column of x (i.e. from mcols(x)$qualities). If x has no qualities metadata column and the qualities argument is omitted, then the fake quality ';' is assigned to each letter in x and written to the FASTQ file. If x is a QualityScaledXStringSet and qualities is not defined, the qualities contained in x are used automatically.

    14. objname	:- The name of the serialized object.

    15. dirpath	:- The path to the directory where to save the serialized object.

    16. save.dups	:- TRUE or FALSE. If TRUE then the Dups object describing how duplicated elements in x are related to each other is saved too. For advanced users only.

    17. verbose	:- TRUE or FALSE.



## Aligning the input sequences 

AlignSeqs() function performs profile-to-profile alignment of multiple unaligned sequences following a guide tree.

```{r}
#AlignSeqs(myXStringSet,
         #guideTree = NULL,
         #iterations = 2,
         #refinements = 1,
         #gapOpening = c(-18, -14),
         #gapExtension = c(-3, -2),
         #useStructures = TRUE,
         #structures = NULL,
         #FUN = AdjustAlignment,
         #levels = c(0.9, 0.7, 0.7, 0.4, 10, 5, 5, 2),
         #alphabet = AA_REDUCED[[1]],
         #processors = 1,
         #verbose = TRUE,
         #...)

```


         
Arguments -

    1. myXStringSet	:- An AAStringSet, DNAStringSet, or RNAStringSet object of unaligned sequences.

    2. guideTree :-	Either NULL or a dendrogram giving the ordered tree structure in which to align profiles. If NULL then a guide tree will be automatically constructed based on the order of shared k-mers.

    3. iterations	:- Number of iteration steps to perform. During each iteration step the guide tree is regenerated based on the alignment and the sequences are realigned.

    4. refinements :-	Number of refinement steps to perform. During each refinement step groups of sequences are realigned to rest of the sequences, and the best of these two alignments (before and after realignment) is kept.

    5. gapOpening	:- Single numeric giving the cost for opening a gap in the alignment, or two numbers giving the minimum and maximum costs. In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.

    6. gapExtension :- Single numeric giving the cost for extending an open gap in the alignment, or two numbers giving the minimum and maximum costs. In the latter case the cost will be varied depending upon whether the groups of sequences being aligned are nearly identical or maximally distant.

    7. useStructures :-	Logical indicating whether to use secondary structure predictions during alignment. If TRUE (the default), secondary structure probabilities will be automatically calculated for amino acid and RNA sequences if they are not provided (i.e., when structures is NULL).

    8. structures	:- Either a list of secondary structure probabilities matching the structureMatrix, such as that output by PredictHEC or PredictDBN, or NULL to generate the structures automatically. Only applicable if myXStringSet is an AAStringSet or RNAStringSet.

    9. FUN :-	A function to be applied after each profile-to-profile alignment.  
    
    10. levels	:- Numeric with eight elements specifying the levels at which to trigger events. 

    11. alphabet	:- Character vector of amino acid groupings used to reduce the 20 standard amino acids into smaller groups. Alphabet reduction helps to find more distant homologies between sequences. A non-reduced amino acid alphabet can be used by setting alphabet equal to AA_STANDARD. Only applicable if myXStringSet is an AAStringSet.

    12. processors	:- The number of processors to use, or NULL to automatically detect and use all available processors.

    13. verbose :- Logical indicating whether to display progress.

    14. ... :- Further arguments to be passed directly to AlignProfiles, including perfectMatch, misMatch, gapPower, terminalGap, restrict, anchor, normPower, standardize, substitutionMatrix, and structureMatrix.



## Viewing aligned Sequences in a Web Browser

BrowseSeqs() function Opens an html file in a web browser to show the sequences in an XStringSet

```{r}
#BrowseSeqs(myXStringSet,
           #htmlFile = tempfile(fileext=".html"),
           #openURL = interactive(),
           #colorPatterns = TRUE,
           #highlight = NA,
           #patterns = c("-", alphabet(myXStringSet, baseOnly=TRUE)),
           #colors = substring(rainbow(length(patterns),
                              #v=0.8, start=0.9, end=0.7), 1, 7),
           #colWidth = Inf,
           #title = "",
           #...)
           
```


Arguments -

    1. myXStringSet	:- A XStringSet object of sequences.

    2. htmlFile	:- Character string giving the location where the html file should be written.

    3. openURL :-	Logical indicating whether the htmlFile should be opened in a web browser.

    4. colorPatterns :- Logical specifying whether to color matched patterns, or an integer vector providing pairs of start and stop boundaries for coloring.

    5. highlight :-	Numeric specifying which sequence in the set to use for comparison or NA to color all sequences (default). If highlight is 0 then positions differing from the consensus sequence are highlighted.

    6. patterns	:- Either an AAStringSet, DNAStringSet, or RNAStringSet object, a character vector containing regular expressions, a list of numeric matrices, or NULL. (See details section below.)

    7. colors	:- Character vector providing the color for each of the matched patterns. Typically a character vector with elements of 7 characters: “#” followed by the red, blue, green values in hexadecimal (after rescaling to 0 ... 255). Ignored when patterns is a list of matrices.

    8. colWidth	:- Integer giving the maximum number of nucleotides wide the display can be before starting a new page. Must be a multiple of 20 (e.g., 100), or Inf (the default) to display all the sequences in one set of rows.

    9. title :-	Character string denoting a title that should appear at the top of the output or "" (the default) for no title.

    10. ...	:- Additional arguments to adjust the appearance of the consensus sequence at the base of the display. Passed directly to ConsensusSequence for an AAStringSet, DNAStringSet, or RNAStringSet, or to consensusString for a BStringSet.  
    
    
    
## Reading aligned sequence files in  fasta format  

read.phyDat() function reads a file in fasta format. This format are used to store nucleotide or protein multiple alignments.

```{r}
# alignment <- read.phyDat("seqs_ali.fasta", format = "fasta")

```

Arguments -

    1. file	:- A file name specified by either a variable of mode character, or a double-quoted string.

    2. format	:- File format of the sequence alignment (see details). Several popular formats are supported: "phylip", "interleaved", "sequential", "clustal", "fasta" or "nexus", or any unambiguous abbreviation of these.

    3. type	:- Type of sequences ("DNA", "AA", "CODON" or "USER").

    4. ... :- Further arguments passed to or from other methods.

    5. x :-	An object of class phyDat.

    6. colsep	:- A character used to separate the columns (a single space by default).

    7. nbcol :- A numeric specifying the number of columns per row (-1 by default); may be negative implying that the nucleotides are printed on a single line.



## Generating a distance matrix for the aligned sequences

dist.ml() function computes pairwise distances for an object of class phyDat. dist.ml uses DNA / AA sequences to compute distances under different substitution models.

```{r}
# dist_matrix <- dist.ml(alignment)

```

Arguments -

    1. x :- An object of class phyDat

    2. ratio :- Compute uncorrected ('p') distance or character difference.

    3. exclude :-	One of "none", "all", "pairwise" indicating whether to delete the sites with missing data (or ambiguous states). The default is handle missing data as in pml.

    4. model :-	One of "JC69", "F81" or one of 17 amino acid models see details.

    5. bf	:- A vector of base frequencies.

    6. Q :-	A vector containing the lower triangular part of the rate matrix.
 
    7. k :- Number of intervals of the discrete gamma distribution.

    8. shape :- Shape parameter of the gamma distribution.

    9. ... :- Further arguments passed to or from other methods.
    
    
    
## Creating a Neighbor-Joining tree with edge lengths
   
NJ() function performs the neighbor-joining tree estimation of Saitou and Nei (1987). UNJ is the unweighted version from Gascuel (1997).

```{r}
# nj_tree <- NJ(dist_matrix)


```


Arguments -

    1. x	:- A distance matrix



## Fitting the NJ tree to the alignment using Maximum Likelihood to get edge lengths

pml() function computes the likelihood of a phylogenetic tree given a sequence alignment and a model. optim.pml optimizes the different model parameters. 


```{r}
# fit_nj <- pml(nj_tree, data = alignment)


```

Arguments -
    1. tree	:- A phylogenetic tree, object of class phylo

    2. data	:- An alignment, object of class phyDat
    
  
    
##  Optimizing the NJ tree to refine the edge lengths


optim.pml() function optim.pml optimizes the different model parameters.


```{r}
# fit_nj <- optim.pml(fit_nj, model = "GTR", optGamma = TRUE, optInv = TRUE)


```

Arguments - 

    1. object	:- An object of class pml.

    2. optNni	:- Logical value indicating whether topology gets optimized (NNI).

    3. optBf :- Logical value indicating whether base frequencies gets optimized.

    4. optQ	:- Logical value indicating whether rate matrix gets optimized.

    5. optInv	:- Logical value indicating whether proportion of variable size gets optimized.

    6. optGamma	:- Logical value indicating whether gamma rate parameter gets optimized.

    7. optEdge :- Logical value indicating the edge lengths gets optimized.

    8. optRate :-	Logical value indicating the overall rate gets optimized.

    9. optRooted :- Logical value indicating if the edge lengths of a rooted tree get optimized.

    10. control :- A list of parameters for controlling the fitting process.

    11. rearrangement :- Type of tree tree rearrangements to perform, one of "none", "NNI", "stochastic" or "ratchet"

    12. subs :- A (integer) vector same length as Q to specify the optimization of Q

    13. ratchet.par	:- Search parameter for stochastic search



## Performing bootstrap resampling (e.g., 100 times)

bootstrap.pml() function performs (non-parametric) bootstrap analysis and bootstrap.phyDat produces a list of bootstrapped data sets. plotBS plots a phylogenetic tree with the bootstrap values assigned to the (internal) edges.

```{r}
# bs <- bootstrap.pml(fit_nj, bs = 100, optNni = TRUE, model = "GTR")

```

Arguments -

    1. x :- An object of class pml or phyDat.

    2. bs	:- Number of bootstrap samples.

    3. trees :-	Return trees only (default) or whole pml objects.

    4. multicore :- Logical, whether models should estimated in parallel.

    5. mc.cores	:- The number of cores to use during bootstrap. Only supported on UNIX-alike systems.

    6. tip.dates :-	A named vector of sampling times associated to the tips/sequences. Leave empty if not estimating tip dated phylogenies.

    7. ... :-	Further parameters used by optim.pml or plot.phylo.

    8. FUN :-	The function to estimate the trees.

    9. jumble	:- Logical, jumble the order of the sequences.
    
    
    
## Creating a consensus tree with the bootstrap value
    
```{r}
# consensus_tree <- plotBS(fit_nj$tree, bs, p = 50)  # Use p = 50 to show bootstrap values >= 50%

```

Arguments -

    1. tree	:- The tree on which edges the bootstrap values are plotted.

    2. trees :- A list of trees (object of class "multiPhylo").

    3. type	:- The type of tree to plot, one of "phylogram", "cladogram", "fan", "unrooted", "radial" or "none". If type is "none" the tree is returned with the bootstrap values assigned to the node labels.

    4. method	:- Either "FBP" the classical bootstrap (default), "TBE" (transfer bootstrap) or "MCC" for assigning clade credibilities. In case of "MCC" all trees need to be rooted.

    5. bs.col	:- Color of bootstrap support labels.

    6. bs.adj	:- One or two numeric values specifying the horizontal and vertical justification of the bootstrap labels.

    7. digits	:- Integer indicating the number of decimal places.

    8. p :- Only plot support values higher than this percentage number (default is 0).

    9. frame :-	A character string specifying the kind of frame to be printed around the bootstrap values. This must be one of "none" (the default), "rect" or "circle".

    10. tol	:- A numeric value giving the tolerance to consider a branch length significantly greater than zero.

    11. sep	:- Seperator between the different methods.

    12. ...	:- Further parameters used by plot.phylo.



## Plotting the tree with bootstrap values

ggtree() function is to draw phylogenetic tree from phylo object

```{r}
# library(ggtree)

# To get node numbers for highlighting particular node of organism of interest/query 

# p <- ggtree(consensus_tree, cex = 0.1, branch.length = 0.1, layout = "rectangular") +
  #geom_tiplab(hjust = -0.00001, vjust = 0.55) +  # Add tip labels
  #xlim(c(0, 0.4)) +  # Adjust x-axis to control branch lengths
  #ggtitle('Phylogenetic tree of IRLC legumes with bootstrap values')


#p1 <- ggtree(consensus_tree, cex = 0.1, branch.length = 0.1, layout = "rectangular") +
  #geom_tiplab(hjust = -0.00001, vjust = 0.55) +  # Add tip labels
  #xlim(c(0, 0.4)) +  # Adjust x-axis to control branch lengths
  #ggtitle('Phylogenetic tree of IRLC legumes with bootstrap values') +
  #geom_nodelab(size = 2.7, hjust = 1.65, vjust = -0.5)+ # Display the bootstrap values on the tree
  #geom_highlight(node = 18, fill = "green", alpha = 0.4,extend = 0.075) # highlight particular node


```


Arguments -

    1. tr	:- phylo object

    2. mapping :- Aesthetic mapping

    3. layout	:- One of 'rectangular', 'dendrogram', 'slanted', 'ellipse', 'roundrect', 'fan', 'circular', 'inward_circular', 'radial', 'equal_angle', 'daylight' or 'ape'

    4. open.angle	:- Open angle, only for 'fan' layout

    5. mrsd	:- Most recent sampling date

    6. as.Date :- Logical whether using Date class in time tree

    7. yscale	:- y scale

    8. yscale_mapping :- yscale mapping for category variable

    9. ladderize :-	Logical (default TRUE). Should the tree be re-organized to have a 'ladder' aspect?

    10. right	:- Logical. If ladderize = TRUE, should the ladder have the smallest clade on the right-hand side? See ape::ladderize() for more information.

    11. branch.length	:- Variable for scaling branch, if 'none' draw cladogram

    12. root.position	:- Position of the root node (default = 0)

    13. xlim :-	x limits, only works for 'inward_circular' layout

    14. layout.params	:- List, the parameters of layout, when layout is a function.

    15. hang :- Numeric The fraction of the tree plot height by which labels should hang below the rest of the plot. A negative value will cause the labels to hang down from 0. This parameter only work with the 'dendrogram' layout for 'hclust' like class, default is 0.1.

    16. ...	:- Additional parameter some dot arguments: nsplit integer, the number of branch blocks divided when 'continuous' is not "none", default is 200.


