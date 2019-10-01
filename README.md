# GFFCompare
Python script which finds matching CDS of different GFF files. 

Intended usage is to compare two gff files which should be near identical, e.g. two GFF files generated by different pipelines but from the same sequence

## Input
* **Sample GFF File**

* **Template GFF File**

Sample GFF file is compared against the Template GFF one.

## Output

* **Strain Stats File**: Contains summary counts of features and mimatches across sample and template files. 
  
   Filename will be TemplateGFFFileName_StrainStats.txt
   An example output would be:
   ```
   Strain	Num_CDS_Automated	Num_CDS_Deposited	Num_Genes_Automated	Num_Genes_Deposited	Num_mRNA_Automated	Num_mRNA_Deposited	Other/MiscFeatures_Automated	Other/MiscFeatures_Deposited	Total_Features_Automated	Total_Features_Deposited	Perfect_Matches	No_Matches	Total_Mismatches	Mismatches_Discarded	Mismatches_Logged	Total_Mismatch_lenght(Logged_Only)	Avg_Mismatch_Lenght
   AY446894.2.fasta.gff	187	188	169	173	168	20	2	132	526	513	183	0	4	3	1	48	48
   ```

* **Mismatch File**: Logs detected mismatches. Mismatch treshold is 3 by default, can be changed via the -m arguement.

   Filename will be TemplateGFFFileName_Mismatch.txt
   An example output would be:
   ```
   Logging mismatches bigger than 3

   Mismatch Origin	ID=UL119_cds2	168891	169263	
   Possible Match	ID=cds-AAR31662.1	168891	169311	
   Mismatch Length 	Start: 0	Stop: 48
   ```

* **Protein/Gene Counts File**:

   Filename will be TemplateGFFFileName_GeneCounts.txt
   A truncated example output would be:
   ```
   Gene	Num_CDS_Automated	Num_CDS_Deposited	Num_Genes_Automated	Num_Genes_Deposited	Num_mRNA_Automated	Num_mRNA_Deposited	Other/MiscFeatures_Automated	Other/MiscFeatures_Deposited	Logged_Mismatches(Sample)
   IRS1	1	1	1	1	1	0	0	0	0
   RL10	1	1	1	1	1	0	0	0	0
   ...
   ```

* **Log File**: Contains general logs of gffcompare's process. Can be very long. Useful for debugging.

   Filename will be TemplateGFFFileName_Log.txt
   A truncated example output would be:
   ```
   Reading ./outputFolder/AY446894.2.fasta_annotation.gff as sample file... 
   GGF Sample file read. 
   Total lines: 526
   CDC lines: 187
   mRNA lines: 168
   Gene lines: 169
   Misc. Features lines: 2
   Other lines: 0
 
   Reading ./manual/AY446894.2.fasta.gff as template file... 
   GFF Template File read.
   Total lines: 513
   CDC lines: 188
   mRNA lines: 20
   Gene lines: 173
   Misc. Features lines: 0
   Other lines: 132
 
   Beginning comparison of CDS features...

   Searching for matching Sample CDS: ('AY446894.2', 'exonerate', 'CDS', '165858', '166799', '.', '-', '.', 'ID=UL116_cds1;Parent=UL116_mRNA;Name=UL116.1;Product=UL116\n')
   Perfect match with Template CDS: ('AY446894.2', 'Genbank', 'CDS', '165858', '166799', '.', '-', '0', 'ID=cds-AAR31660.1;Parent=gene-HHV5wtgp102;Dbxref=NCBI_GP:AAR31660.1;Name=AAR31660.1;gbkey=CDS;gene=UL116;locus_tag=HHV5wtgp102;product=protein UL116;protein_id=AAR31660.1\r\n')
   Start positions: 165858 165858
   Stop positions: 166799 166799
   ...
   ```

## Arguments
gffcompare.py S T [-i -m -p]

### Required

* S: Path of sample GFF file, e.g. ~/directory/SampleGFF.gff

* T: Path of template GFF file,  e.g. ~/directory/TemplateGFF.gff
                    
### Optional

* -m    
  **Mismatch Permissibility**: The max length of mismatches that will NOT be logged in mismatch file. Default is 3.
                    
* -i                 
  **Iteration Number**: Max times to run fuzzy filter. Each iteration increases range by mismatch permissibility (default 3). Default is 10

* -p    
  **Gene Name List File**: Directory of file with gene names (1 per line) to be used in generating the counts table. Default is path/proteinList.txt. A Human Cytomegalovirus Protein List is provided as an example.
   
   
## TODO (Low Priority)

Rename headers in Gene Counts file (as well as name of file itself) for consistency

Short blurb on how matching works

Diagram???
   
