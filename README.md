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
* GATK CreateSequenceDictionary 

## Diagram
<p align="center"> 
<img src="https://user-images.githubusercontent.com/73086054/189310509-375fea4f-11fb-41ca-ba52-90760e9a5aa3.png" width="80%">
</p> 

## User guide
**1. Set up**

Clone this repository by running:
```
git clone https://github.com/Sydney-Informatics-Hub/IndexReferenceFasta-nf.git
cd IndexReferenceFasta-nf
``` 

**2. Generate indexes**  

Users can specify which index files to create by using the `--bwa`, and/or `--gatk` flags. GATK and BWA indexes are optional, while Samtools is run by default. Run the pipeline with:

```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk 
```

If you are running the pipeline on NCI Gadi or Pawsey's Nimbus cloud you should specify this with the `-profile` flag at runtime. This will allow you to use Singularity to run containers at Gadi and Docker to run the containers at Nimbus. 

To run the pipeline at NCI Gadi:

```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk -profile gadi --whoami <us1111> --pbs_account <aa00>
``` 

To run the pipeline at Pawsey's Nimbus cloud:

```
nextflow run main.nf --ref /path/to/ref.fasta --bwa --gatk -profile nimbus
```

Infrasturcture-specific config files can be found in `config/`

## Benchmarking

### Human hg38 reference assembly @ Pawsey's Nimbus (NCPU/task = 1)
|task_id|hash     |native_id|name          |status   |exit|submit |duration  |realtime  |%cpu   |peak_rss|peak_vmem|rchar  |wchar  |
|-------|---------|---------|--------------|---------|----|-------|----------|----------|-------|--------|---------|-------|-------|
|3      |27/33fffc|131621   |samtools_index|COMPLETED|0   |55:44.9|12.2s     |12s       |99.20% |6.3 MB  |11.8 MB  |3 GB   |19.1 KB|
|1      |80/f03e46|131999   |gatk_index    |COMPLETED|0   |55:46.7|22.6s     |22.3s     |231.90%|3.8 GB  |37.1 GB  |3.1 GB |726 KB |
|2      |ea/e29535|131594   |bwa_index     |COMPLETED|0   |55:44.9|1h 50m 16s|1h 50m 15s|99.50% |4.5 GB  |4.5 GB   |12.1 GB|8.2 GB |

## Workflow summaries

### Metadata
|metadata field     | workflow_name / workflow_version  |
|-------------------|:---------------------------------:|
|Version            | 1.0                               |
|Maturity           | stable		                |
|Creators           | Georgie Samaha                    |
|Source             | NA                                |
|License            | GPL-3.0 license                   |
|Workflow manager   | NextFlow                          |
|Container          | None                              |
|Install method     | Manual                            |
|GitHub             | Sydney-Informatics-Hub/IndexReferenceFasta-nf                                |
|bio.tools          | NA                                |
|BioContainers      | NA                                | 
|bioconda           | NA                                |

### Component tools

* samtools/1.15.1
* gatk/4.3.0.0 
* bwa/0.7.17

### Required (minimum) inputs/parameters

* A reference genome file in fasta format.

## Additional notes

### Help/FAQ/Troubleshooting

* A subset fasta file for testing is available in `testData/` 

## Acknowledgements/citations/credits
### Authors 
- Georgie Samaha (Sydney Informatics Hub, University of Sydney)   

### Acknowledgements 

- This pipeline was built using the [Nextflow DSL2 template](https://github.com/Sydney-Informatics-Hub/Nextflow_DSL2_template).  
- Documentation was created following the [Australian BioCommons documentation guidelines](https://github.com/AustralianBioCommons/doc_guidelines).  

### Cite us to support us! 
Acknowledgements (and co-authorship, where appropriate) are an important way for us to demonstrate the value we bring to your research. Your research outcomes are vital for ongoing funding of the Sydney Informatics Hub and national compute facilities. We suggest including the following acknowledgement in any publications that follow from this work:  

The authors acknowledge the technical assistance provided by the Sydney Informatics Hub, a Core Research Facility of the University of Sydney and the Australian BioCommons which is enabled by NCRIS via Bioplatforms Australia. 
