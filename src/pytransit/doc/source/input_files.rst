.. _input_files:

File Formats for Transit
=======================

.. _combined_wig_link:

Combined_wig Files
------------------

Transit now supports a new file format called 'combined_wig' which basically
combines multiple wig files into one file (with multiple columns).  This is
especially useful for managing larger collections of datasets.

Combined_wig files can created through the Transit GUI
(in the menu: **Pre-processing -> Export -> Combined_wig**), or via the command line:

::

  Usage:

  > python3 transit.py export combined_wig --help

  usage: python3 transit.py export combined_wig <comma-separated .wig files> <annotation .prot_table> <output file> [-n <norm>]

  Example:

  > python3 transit.py export combined_wig Rv_1_H37RvRef.wig,Rv_2_H37RvRef.wig,Rv_3_H37RvRef.wig H37Rv.prot_table clinicals_combined_TTR.cwig

You can specify the normalization method you want to use with a flag.
TTR is the default, but other relevant normalization options would be 'nonorm'
(i.e. preserve raw counts) and 'betageom' (this corrects for skew, but is slow).

The format of a combined_wig is simply a multi-column file with
the first column being the coordinates of TA sites, followed by 
N columns of counts (for N samples), possibly with a final column indicating
the gene annotation information.
A combined_wig file can have header lines, prefixed by '#'.

Importantly, a combined_wig file must include sample identifiers
(filenames) that are prefixed with **"#File: "**.  These header lines
are automatically included by the 'export combined_wig' Transit
command the creates a combined_wig from multiple .wig files.  These "File:"
header lines should be given in the same order as the sample columns,
and the filenames are used as identifiers to cross-reference
information about each sample in the metadata file (see below).

::

 #command: python3 ../src/transit.py export combined_wig /Users/example_data/glycerol_rep1.wig,/Users/example_data/glycerol_rep2.wig,/Users/example_data/cholesterol_rep1.wig,/Users/example_data/cholesterol_rep2.wig,/Users/example_data/cholesterol_rep3.wig H37Rv.prot_table temp.cwig
 #genome: H37Rv
 #File: src/pytransit/data/glycerol_rep1.wig
 #File: src/pytransit/data/glycerol_rep2.wig
 #File: src/pytransit/data/cholesterol_rep1.wig
 #File: src/pytransit/data/cholesterol_rep2.wig
 #File: src/pytransit/data/cholesterol_rep3.wig
 #normalization: TTR
 #Normalization factors: 1.89 2.01 0.97 1.15 1.06
 #TA  glycerol_rep1  glycerol_rep2  cholesterol_rep1  cholesterol_rep2  cholesterol_rep3  annot
 60	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 72	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 102	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 188	0.0	0.0	0.0	0.0	0.0	Rv0001 (dnaA)
 ...

(DnaA is essential, which is why there are no insertion counts in this example)

Note that the ORF id and gene names are appended for TA sites in CDS regions, for convenience.
If you open a combined_wig file in Excel (with tab as separator character), you will
see the last line of the header (starting with "TAcoord...") provides the input wig filenames as column headers.

The '#genome' is extracted from information in the individual .wig files,
which records the reference genome sequence that was used with TPP; 
the coordinates of TA sites are determined from the genome sequence.
The **reference genome must be the same for all wigs** being combined.


|


.. _metadata_files:

Samples Metadata File
---------------------

The **metadata** file describes the sample ids, filenames,
and conditions they represent (e.g. different media, growth
conditions, knockout strains, animal passaging, etc., or whatever
treatments and controls your TnSeq experiment involves) for all the
samples in a combined_wig file.  

Format of the *samples_metadata* file: a tab-separated file (which you
can edit in Excel) with 3 columns: Id, Condition, and Filename (it
must have these headers). The Condition column should be as specific
as possible, indicating **how to group replicates**.
You can include other columns of info, but
do not include additional rows.  Individual rows can be commented out
by prefixing them with a '#'.  Here is an example of a samples
metadata file: The filenames should match what is shown in the header
of the combined_wig (including pathnames, if present).

Note: the Condition column should have a unique label for each
distinct condition (the same label shared only among replicates).  If
there are attributes that distinguish the conditions (such as strain,
treatment, etc), they could be included as **additional columns**
(e.g. covariates, like Carbon_source or Drug or Days or Strain...).

::

  Id      Condition    Filename
  glyc1   glycerol     /Users/example_data/glycerol_rep1.wig
  glyc2   glycerol     /Users/example_data/glycerol_rep2.wig
  chol1   cholesterol  /Users/example_data/cholesterol_rep1.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep2.wig
  chol2   cholesterol  /Users/example_data/cholesterol_rep3.wig


|

.. _annotation_files:

Genome Annotations (.prot_tables and .gff files)
------------------

The annotation of a genome contains information about genes, such as
coordinates, strand, locus tag, gene name, and functional description.
Transit uses a custom format for annotations called "prot_table"s,
e.g. H37Rv.prot_table.  Prot_tables are **tab-separated text files**
containing the gene information in 9 specific columns:

**Prot_table file format:**

1. gene function description
2. start coordinate
3. end coordinate
4. strand
5. length of protein product (in amino acids)
6. don't care
7. don't care
8. gene name (like "dnaA")
9. ORF id (like Rv0001)

Examples of prot_tables for commonly used genomes can be found at:
`https://orca1.tamu.edu/essentiality/transit/genomes/
<https://orca1.tamu.edu/essentiality/transit/genomes/>`_.

Here is an example (H37Rv.prot_table):
  
::

  chromosomal replication initiation protein 	1	1524	+	507	15607143	885041	dnaA	Rv0001
  DNA polymerase III subunit beta 	2052	3260	+	402	15607144	887092	dnaN	Rv0002
  recombination protein F 	3280	4437	+	385	15607145	887089	recF	Rv0003
  hypothetical protein Rv0004 	4434	4997	+	187	15607146	887088	-	Rv0004
  DNA gyrase subunit B 	5123	7267	+	714	15607147	887081	gyrB	Rv0005
  DNA gyrase subunit A 	7302	9818	+	838	15607148	887105	gyrA	Rv0006
  ... 

  (full file has ~4000 lines)


.. NOTE::

  *It is crucial* that the annotation file (.prot_table) used for
  analyses in Transit corresponds to exactly the same genome sequence
  (.fasta or .fna) that was used to generate the .wig files with TPP,
  because it is used to determine which TA sites are contained in which
  genes (by coordinates). For example, **H37Rv.fna** is paired with
  **H37Rv.prot_table**, both derived from GenBank sequence NC_000962.2.

In many cases, users might often obtain annotations for their genome
in **.gff (or .gff3)** file format, such as downloaded from NCBI.  .gff
files contains essentially the same information about genes.  However,
there is a bit more flexibility in the .gff file format (especially in
the tags used in the right-most column), and the information about
genes is not always encoded in a uniform way, making it difficult to
use arbitrary .gff filess for analyses in Transit.  
Therefore, there is a
simple procedure in Transit to convert a .gff file to .prot_table
format (via GUI or command-line).  This
step only has to be done once, and then the .prot_table can be used
for all subsequent analyses in Transit.
(The routine specifically looks for the 'locus_tag', 'gene', and 'product'
tags in info field of CDS records.)

::

  > python3 transit.py convert gff_to_prot_table <input.gff_file> <output.prot_table>




|

.. rst-class:: transit_sectionend
----
