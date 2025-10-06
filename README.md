IndexReferenceFasta-nf
===========

  - [Description](#description)
  - [Diagram](#diagram)
  - [User guide](#user-guide)
  - [Benchmarking](#benchmarking)
  - [Workflow summaries](#workflow-summaries)
      - [Metadata](#metadata)
      - [Component tools](#component-tools)
      - [Required (minimum)
        inputs/parameters](#required-minimum-inputsparameters)
  - [Additional notes](#additional-notes)
  - [Help/FAQ/Troubleshooting](#helpfaqtroubleshooting)
  - [Acknowledgements/citations/credits](#acknowledgementscitationscredits)

---

## Description
This is a flexible pipeline for generating common reference genome index files for WGS data analysis. IndexReferenceFasta-nf is a Nextflow (DSL2) pipeline that runs the following tools using either Docker or Singularity to run containerised software for:
* Samtools faidx
* BWA index
* BWA-mem2 index
* GATK CreateSequenceDictionary 

## User guide
**1. Set up**

Clone this repository by running:
```
git clone https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf.git
cd IndexReferenceFasta-nf
``` 

**2. Generate indexes**  

Users can specify which index files to create by using the `--bwa`, `--bwamem2`, `-samtools`, and/or `--gatk` flags. Run the pipeline with:

```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk
```

If you are running the pipeline on NCI Gadi you should specify this with the `-profile` flag at runtime. This will allow you to use Singularity to run containers at Gadi.

**Standard**

To run the pipeline on your own system, you will need to have Nextflow installed. You can adjust the `standard.config`  configuration file depending on your own system needs. Currently it runs containers with Singularity. You can test your customised config file using the test fasta available in `testData`.

To run the pipeline with the `standard.config`, run the following: 

```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk -profile standard 
```

**NCI Gadi HPC**  

To run the pipeline at NCI Gadi, first load the Gadi-specific Nextflow installation:

```
module load nextflow
```

Then run the pipeline:
```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk -profile gadi --whoami <us1111> --pbs_account <aa00>
``` 

### Metadata
|metadata field     | workflow_name / workflow_version  |
|-------------------|:---------------------------------:|
|Version            | 2.0                               |
|Maturity           | stable		                |
|Creators           | Georgie Samaha                    |
|Source             | NA                                |
|License            | GPL-3.0 license                   |
|Workflow manager   | NextFlow                          |
|Container          | None                              |
|Install method     | NA                            |
|GitHub             | Sydney-Informatics-Hub/IndexReferenceFasta-nf                                |
|bio.tools          | NA                                |
|BioContainers      | NA                                | 
|bioconda           | NA                                |

### Component tools

* samtools/1.15.1
* gatk/4.3.0.0 
* bwa/0.7.17
* bwa-mem2/2.3

### Required (minimum) inputs/parameters

* A reference genome file in fasta format.

## Additional notes

### Help/FAQ/Troubleshooting

* A subset fasta file for testing is available in [`testData/`](https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf/tree/main/testData) 

## Acknowledgements/citations/credits
### Authors 
- Georgie Samaha (Sydney Informatics Hub, University of Sydney)   

### Acknowledgements 

- This pipeline was built using the [Nextflow DSL2 template](https://github.com/Sydney-Informatics-Hub/Nextflow_DSL2_template).  
- Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).  

### Cite us to support us! 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:  

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia. 
