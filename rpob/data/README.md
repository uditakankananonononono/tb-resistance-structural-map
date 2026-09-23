# README

## Summary

This document contains a history of how, why and when the CRyPTIC final data tables were updated.

For more information on each table and how they are related, please see DATA_SCHEMA.pdf

Philip W Fowler, 14 Oct 2018

## Unique keys

Most tables rely on UNIQUEID, perhaps in combination with one or more other fields (such as READINGDAY or DRUG), to make a unique index. The unique index (i.e. primary key) of each data table is marked on the DATA_SCHEMA with a key symbol and the fields are drawn in /italics/.

This is an example of a UNIQUEID and is closely related to the UniqueID that is passed down from Clockwork in the VCF metadata.

`site.03.subj.T334.lab.T334.iso.1`

Its format is the

`site.SITEID.subj.SUBJID.lab.LABID.iso.ISOLATENO`

SITEID is the "%02d" two-digit string identifying the site e.g. 05 is Peru.
SUBJID is a string identifying the patient/source
LABID is a string identifying the sample taken from that SUBJID. For several labs this is the same as SUBJID
ISOLATENO is a integer describing which isolate from the sample was used. This is almost always 1.

One way of thinking about UNIQUEID is therefore the cultured sample in the pipette before it is (a) inoculated onto a plate and/or (b) sent for sequencing.

PHENOTYPE data therefore requires in addition READINGDAY to identify a unique reading of the plate, whilst GENOTYPE data therefore requires SEQ_REPS, which is the 'seqencing repeat' and is usually 1, but because if there is more than one repeat, Clockwork simply takes all the short reads it can find, it is stored as a string e.g. "1_2_3" means there were 3 sequencing repeats and Clockwork has used data from all of them. If there is more than one VCF per UNIQUEID, we therefore take the one with the longest SEQ_REPS.

Most tables contain SITEID as well as UNIQUEID since although the former can be easily derived from the latter, since we often look at by site statistics, it is there to just save time.

## Documentation

You'll find in this folder several important documents

```
DATA_SCHEMA       This shows the relationships between the key tables and their primary keys. If you load a .pkl version it will be indexed accordingly.
DATA_NUMBERS      A simple high-level view of the *current* number of samples represented in the tables.
DATA_WALKTHROUGH  Takes you through all the tables and explains what they contain and how they can be used
```

## WGS_PREDICTION_STRING

In `GENOMES` table; noted here as is useful! The drug order is

```
["RIF","INH","PZA","EMB","AMI","KAN","LEV","MXF","ETH","PAS","RFB","LZD","BDQ","DLM","CFZ"]
```

## History

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
