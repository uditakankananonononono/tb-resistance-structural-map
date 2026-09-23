
# Release Notes

## cryptic-tables-v3.4.0, 21 May 2025

Main changes

* Redownloaded all CRyPTIC 96-well plate data from clires2.org. This picked up an additional 425 samples from National University of Singapore
* Replaced AMyGDA with TMAS -- this is a convolutional neural network that has been trained on a set of CRyPTIC images to read MICs from the photographs taken by Vizion. See https://doi.org/10.1101/2025.02.14.638231. This has increased the proportion of MICs in which we have a high confidence (because at least two of the three measurements methods agree) from 79.2% to 88.7%.
* Added two tables with the 4.75 million classifications made by the BashTheBug volunteers of the MIC for each drug on all the CRyPTIC/UKMYC 96-well plates.

Future work

* process the FASTQ files from NUS so the new phenotypes have matching genetics

## cryptic-tables-v3.3.0, 24 April 2025

Main changes

* Added WGS for additional 4,832 samples that were processed through Pathogena / sp3dev in early April 2025

Future work

* revamp the PHENOTYPE_QUALITY for MICs measured via the UKMYC plates using the new TMAS machine learning model. Expect this will replace AMyGDA and increase the proportion of MICs we can annotate as being High confidence due to agreement between at least two independent measurement methods. This is ongoing.

## cryptic-tables-v3.2.0, 8 April 2025

In this version 

* Some additional pDST data has been included that Kerri Malone identified as part of the `cryptic-catalogues` project that Dylan Adlard has taken over. An additional 9,160 samples with 54,910 pDST measurements have been added to DST_MEASUREMENTS. Most were already identified in WGS_SAMPLES but a further 134 samples have been added (labelled Group 6)
* In addition a bug was identified whereby 1,789 samples with existing pDST measurements that included the ENA sample accession as part of their UNIQUEID had not been picked up and added to WGS_SAMPLES.

Note that the larger `pandas.Dataframes` are only stored as parquet files for now -- you may need to e.g. `pip install pyarrow` to load these using `pandas`.

Future work
* process the 1,923 samples in groups 6 and 7 which have been added in this version. Note that these don't have `COUNTRY_CODE` in `DST_SAMPLES` which may need revisiting.
* process the 2,909 from CRyPTIC-v3.0 which didn't upload, most likely because they have three FASTQ files in the ENA 
* so far we've only processed samples which have pDST and WGS; process some/all samples which have WGS but not pDST (useful for e.g. building a very large tree) 
* add more documentation and a schema
* provide CSV files
* begin to expand to include NTMs (a lookup for the NTMs held in the Oxford tenancy is now included -- these will be become public upon publication of the Mycobacterial competitive manuscript) 
* deal with ENA samples that have multiple run accessions; these have been skipped over for now
* begin to add ONT processed samples (this will require some refactoring to e.g. account for one sample being sequenced two ways)


## cryptic-tables-v3.1.0, 28 Feb 2025

This is referred to as "release three" as all these samples have had their FASTQ files downloaded and processed through EIT Pathogena so have been treated fundamentally differently to the preceeding datasets.

This minor version bump combines two datasets; the first is the set of samples uploaded to Pathogena over Christmas 2024. Analysis showed there were several problems with this dataset. First as CRyPTIC UNIQUEID had been used to label the FASTQ files and some UNIQUEIDs had multiple ENA run accessions, there was no way of knowing which ENA run accession had been used. Second, some samples were uploaded but failed to reach a final complete status and a third group were not uploaded to begin with. These were labelled as groups 1,2 and 3, respectively. To these were added two further groups of samples by comparing to work Kerri Malone had carried out in identifying a validation dataset for building a resistance catalogue. These were groups 4  (a small number of samples Dylan Adlard identified as being in the training set which hadn't yet been processed) and 5 (a larger group of samples that consitute a genuine validation dataset.)

These two datasets have been combined for this release increasing the number of samples from 42,637 to 49,263 with both pDST and WGS data. We've introduced a new table `WGS_SAMPLES` (renamed from `SAMPLES` to avoid confusion with `DST_SAMPLES`) which is a subset of `ENA_LOOKUP` and describes which ENA run accessions were downloaded and processed through EIT Pathogena. In here is a column `DATASET` which describes which dataset the sample belongs to. Note that this tables contains samples that we tried to process and hence not all worked due to a variety of reasons. The 25 samples which don't match anything I don't know where they come from so decided to leave as unknown. The values are

```
CRyPTIC-v1.0    37984
CRyPTIC-v2.0     5737
CRyPTIC-v3.0     8896
NaN                25
```
where CRyPTICv1.0 means the sample was in CRyPTIC Release One which was handed to FIND/Seq&Treat to build WHOv1 and also the few samples Dylan Adlard identified as missing (group 4). CRyPTICv2.0 are additional CRyPTIC samples that gained WGS, pDST or both after the data freeze in April 2020. This includes about 1100 samples from NICD that are enriched for BDQ resistance.  Finally CRyPTICv3.0 is the set of 8,896 samples which Kerri Malone identified as being genuine validation samples (group 5), however as shown below only 5,827 reached complete. 

```
status	        cannot      cannot
            	assemble    speciate	complete	not uploaded
dataset				
CRyPTIC-v1.0	    9	        14	      37893	          68
CRyPTIC-v2.0	  129	        83	       5523	           2
CRyPTIC-v3.0	  158	         2	       5827	        2909
NaN	                1	         2	         21	           1
```

## cryptic-tables-v2.1.2, 23 Jul 2024

1. ENA_LOOKUP: thanks to some heavy lifting by @Jeff the ENA_LOOKUP table now contains the mapping between CRyPTIC `UNIQUEID` and ENA sample, run accession numbers. It is potentially one-to-many because (i) a few hundred samples may have been uploaded twice (or the site incorrectly didn’t tell us it was a sequencing repeat, hence have left) or (ii) many of the “ENA” samples that were downloaded from the ENA have multiple run accessions which is reflected in multiple SEQ_REPs in the GENOMES table (this means the FASTQ files were concatenated together before processing by Martin and Jeff). 

The idea is one can now identify samples with specific mutations / phenotype profiles, and then by linking to this table and unpacking the `fastq_ftp` column which has a semicolon delimiter and looks like 

```
ftp.sra.ebi.ac.uk/vol1/fastq/ERR218/006/ERR2184206/ERR2184206_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR218/006/ERR2184206/ERR2184206_2.fastq.gz
```

and after unpacking one can programmatically get the FASTQ files in turn via something like

```
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR218/006/ERR2184206/ERR2184206_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR218/006/ERR2184206/ERR2184206_2.fastq.gz
```

More ambitiously, it will also let us process (nearly) all of CRyPTIC back through GPAS should we wish to do that in the future (to e.g. support a WHO v3 resistance catalogue).

2. Bug in SITE 07 naming: PHE use `/` and `.` in their lab identifiers which breaks code, file paths etc. — in the previous version they were converted to `_` in the phenotypes but not in the genetics tables which meant they wouldn’t join. This has now been fixed so the number of samples where there are genetics and pDST has increased slightly (these were nearly all MGIT first line only samples)

## cryptic-tables-v2.1.1, 20 Feb 2024

Feature: epistasis rules have been added to our version of WHO2 (should improve performance of several drugs) and Jeremy has updated gnomonicus etc so they can be parsed.

This also caused us to reconsider how we report the effect of 

(i) 	minor alleles in resistance genes
(ii) 	all variants detected in mmpL5

Because we are explicitly reporting any evidence of all minor alleles in the resistance genes (by overriding the MIN_FRS filter in Clockwork VCFs) the tables contain a large number of minor alleles, yet none of these are associated with resistance in WHO2. Previously they were reported as having an Unknown phenotype but on reflection it is more appropriate to report these as Susceptible.

Whilst it is categorised as a Tier 1 gene in WHO2, mmpL5 isn’t what we would call a resistance gene since it has no variants associated with resistance. However, we need to detect and report all genetic variation in mmpL5 so we can later apply the new epistasis rules (briefly, any loss of function mutation in mmpL5 overrides any resistance-associated mutation in Rv0678 leading to an overall prediction of Susceptible). Hence all genetic variation in mmpL5 is reported as Susceptible as Unknown only makes sense if there is a reasonable likelihood that the observed mutation could be associated with Resistance.

Our interpretation of the WHO2 catalogue has therefore also changed but this is not yet publicly available.

## cryptic-tables-v2.1.0, 20 Feb 2024

This update focusses on the resistance catalogues. Whilst developing the first draft of the second edition of the WHO catalogue (WHO2) we found a few minor mainly "off by one" bugs in WHO1. Hence all samples have been rescored through both this updated version of WHO1 (`WHO-UCN-GTB-PCI-2021.7` version 1.2) and also the current preliminary version of WHO2 (`WHO-UCN-GTB-PCI-2023.5` version 2.0). Why is it called preliminary? Well the epistasis rules haven't been added yet as these require code changes to `piezo`. As part of this the `WGS_PREDICTION_STRING` in `GENOMES` has been recalculated using WHO2 and obviously the `EFFECTS` and `PREDICTIONS` tables now contain the results for two catalogues so you need to select the catalogue you want to use. Also the `UKMYC_SUBJECTS` and `UKMYC_SAMPLES` tables were previously included in the `DATA_SCHEMA` by mistake so these have been removed.

## cryptic-tables-v2.0.1, 4 Jan 2024

Bug fix - the UNIQUEID in the Genetics tables (VARIANTS, MUTATIONS, EFFECTS and PREDICTIONS but not GENOMES) was incorrect. It had some extra information at the end which prevented joining. This has now been fixed. Nothing else changed.

## cryptic-tables-v2.0.0

Please note the data schema has been simplified with four groupings of tables; these are
1. Genetics. Information on all samples for which we have paired Illumina FASTQ files and have been successfully processed by Clockwork.
2. UKMYC phenotypes. Minimum inhibitory concentrations measured from UKMYC5/6 plates as part of the CRyPTIC project. Not all samples were sequenced and some samples/sites had processing issues. Includes AMyGDA and BashTheBug measurements.
3. All phenotypes. Superset of 2. i.e. includes R/S binary result derived using CRyPTIC ECOFFs as well as any and all phenotypic information we have for any sample. This can be MGIT960, LJ, MYCOTB plate and others. Some samples have multiple SIR measurements for the same drug using different methods.
4. Reference. Series of lookup tables to convert e.g. 3 letter drug code into proper name. Now includes detailed ENA lookup for most (not yet all) samples allowing the original FASTQ files to be automatically downloaded from the ENA.

Detailed notes below. To save disc space, only gzipped pickles of Pandas DataFrames are provided for now.

### Genetics
* all "Genetics" tables now used Clockwork v0.12.4, however not all samples present in v1 have been processed e.g. none of the `site.ENA...` samples are present, including those with historic phenotypic data
* all tables have been created using `gnomonicus` which is new and uses `gumpy` and `piezo` which themselves have been rewritten by Jeremy Westhead. Hence the `MUTATIONS` and `VARIANTS` tables in particular are different with fewer Booleans and some columns renamed.
* one key improvement is we have identified potential minor alleles in genes known to be associated with resistance; `gumpy` effectively has been asked to ignore the `MIN_FRS` filter in the VCF file for these genes and instead report this is a `MINOR_ALLELE`. 
* no samples have been run through Mykrobe so species and lineage information is currently missing.
* note also that all samples had their per-sample VCF processed (which is a change compared to CRyPTIC v1 data tables)

### UKMYC phenotypes

* downloaded all CRyPTIC phenotypic data and images uploaded to CliRes after the "Data Freeze" in 2020; all images have been processed with AMyGDA and all drug lanes have been processed by BashTheBug, hence all MICs have been assigned LOW/MEDIUM/HIGH quality phenotypes where possible

### All phenotypes
* added new CRyPTIC MICs downloaded from CliRes as described above (labelled CRyPTIC2) 
* added any additional DST entered by CRyPTIC labs into CliRes (labelled CLIRES2023)
* added samples with BDQ and LZD resistance from NICD (labelled NICD2022)
* since this table contains phenotypes derived outside of CRyPTIC using broth microdilution plates and these have only been read once, have allowed CRyPTIC MICs with a low quality into this table. Hence all phenotypes now have a `QUALITY`. If `HIGH` this means the measurement has been checked within CRyPTIC by at least two independent methods which agree; if `MEDIUM` it means only one measurement method was used and if `LOW` it means the measurement was checked with CRyPTIC by at least two independent methods, none of which agreed. For the latter the provided MIC is always the one measured by the laboratory scientist.
* corrected some data entry issues in the microtitre data from Seq&Treat in `DST_MEASUREMENTS` e.g. ">0,12" -> ">=0.12"; this allows these to correctly assign a `PHENOTYPE` by joining to `PLATE_LAYOUT`
* due to a mistake pivotting an Excel sheet, there were some rows in `DST_MEASUREMENTS` from Mumbai that had no `MIC` or `PHENOTYPE`; these have been removed as they are potentially confusing

### Reference

* added the `MYCOTB` plate design to `PLATE_LAYOUT` including ECOFFs for any drugs on the UKMYC plates
* removed the `I` category from UKMYC5/6 in `PLATE_LAYOUT` and simply applied the ECOFF as defined in the paper; this means there will be no UKMYC `I` results in `DST_MEASUREMENTS`
* added new `ENA_LOOKUP` table that has ENA study, sample, run accessions for as many samples as possible; includes ENA FTP paths, MD5s and size of FASTQ files in bytes to help with download

## Release One updates

(most recent first)

```
Date            Description
26-Jan-2021     Release 1.1.
03-Jun-2020     Release One. Early access - the tables have been substantially changed! Please have a play and let us know any problems. There are a few glitches to fix which will be done over the next week. Also the ENA and Seq&Treat genetics remain to be added.
06-Apr-2020     BTB finished classifying all (frozen) images; PHENOTYPES and associated tables updated. CLASSIFICATIONS tables updated.
27-Mar-2020     ALL MIC tables. List of UNIQUEIDs in PLATE_READINGS now frozen to help with genetic joint-genotyping by the EBI team.
10-Feb-2020     ALL MIC tables updated to help with identifying data errors. 10,993 matching samples.
28-Jan-2020     PHENOTYPES and PLATE_LAYOUT tables updated to reflect current proposed ECOFFs.
23-Jan-2020     ALL tables updated following complete CliRes rebuild to catch samples renamed to fix some matching errors. More VCFs processed. 10,949 samples matching. ECOFFs slightly altered in PLATE_LAYOUT to reflect current ECOFF paper draft (ATU removed from AMI, KAN, ATU added to MXF, EMB made consistent between plate designs). PREDICTIONS table NOT yet updated with resulting new BINARY_PHENOTYPES (waiting for rescomp to come back online).
12-Dec-2019	ALL MIC tables updated. OTHER_PHENOTYPES modified (LID2015,NEJM2018 rows added, where genetics) 9,094 samples matching. Last update before rebuilding.
18-Nov-2019     ALL tables updated (new 03 vcfs added). 8,885 samples matching.
06-Nov-2019     ALL MIC tables updated, including BASHTHEBUG bug fix making 5k more HIGH MIC measurements available
23-Oct-2019	PHE testing results added to OTHER_PHENOTYPES
15-Oct-2019     All tables updated + few more bugs ironed out. DATA_WALKTHROUGH and DATA_SCHEMA updated
14-Oct-2019     ALL MIC tables updated + few bugs ironed out of tables. SAMPLES and SUBJECTS added.
25-Sep-2019	Huge update; All tables updated and perhaps changed, IR2 per-sample VCFs processed. See email.
24-Jun-2019     All MIC tables updated
12-Jun-2019     All MIC tables updated
27-May-2019	All MIC tables updated
27-Apr-2019     All VCF tables updated with new fishing genes (all mmpS/Ls)
20-Apr-2019     All MIC tables updated.
05-Apr-2019     All MIC tables updated and all ENA VCFs run through and added to genetic tables. CSVs now gzipped.
29-Mar-2019	Added WHO-PHENOTYPES which contain 2,470 samples from the WHO with MGIT/LJ phenotypes, all of which (except 402 which were IonTorrent or both fastqs were not available) will match ENA genomes.
22-Mar-2019     All MIC tables updated.
08-Mar-2019     Large update: All MIC tables updated + added pheno/geno discrepancy boolean flags (e.g. RIF_PGD) + added REF_COVERAGE, ALT_COVERAGE + MINOS_SCORE from Clockwork to MUTATIONS.
10-Feb-2019	Added CATALOGUE_NAME field to the EFFECTS and PREDICTIONS tables. In effect this is the default catalogue, which at the moment is CRYPTICv1.0 = NEJM2018 for RIF,INH,EMB,PZA and LID2015B for other drugs (if available)
08-Feb-2019	All MIC tables updated (for annual report).
04-Feb-2019     All Genetic tables updated using new catalogue which is simply a merge of NEJM2018+LID2015. Will be rapidly adding rows for e.g. ETH, BDQ and re-running but bear with me!
01-Feb-2019     All MIC tables updated.
25-Jan-2019     All MIC tables updated.
18-Jan-2019     All MIC tables updated. MIC measurements from site 13 (Taiwan) added, all match genomes.
14-Jan-2019     All MIC tables updated with new data. Added BASHTHEBUGPRO table.
07-Jan-2019     All tables updated with new data including UKMYC6 plates for the first time. PLATEDESIGN and LOG2MIC columns added.
18-Nov-2018     CliRes-tables updated + recent BashTheBug results added
28-Oct-2018     CliRes-dependent tables updated by PWF to give new totals
18 Oct 2018     Added new columns to PHENOTYPES data giving estimate of quality of reading + update numbers
16 Oct 2018     Added OTHER_PHENOTYPES table with MGIT data sent by site 04 (Mumbai)
16 Oct 2018     BASHTHEBUG tables updated with classifications up to Mon 15 Oct.
14 Oct 2018     cryptic-tables/ created by PWF. Read-only DropBox folder
```


