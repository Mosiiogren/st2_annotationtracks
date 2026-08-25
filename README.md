# StackStorm pack for creating and uploading annotationtrackfiles to Gens

## Installation
```bash
st2 pack install https://github.com/Mosiiogren/Annotationtracks
st2 pack config Annotationtracks
```
### Config file
The following config parameters need to be defined:
- `hosts`: The server the command is running at.
- `cwd`: Path to the folder containing the docker-compose file.
- `outputfileregulatory`: Path to where the regulatory data should be stored.
- `outputfilegene`: Path to where the filtered gene data should be stored.
- `outputfileexon`: Path to where the filtered exon data should be stored.
- `outputfileclusters`: Path to where the clusters should be stored.
- `outputfolderannotationtracks`: Path to where the annotationtrack files should be stored.

Additionally, the following parameter needs to be set up as a st2 key:
```bash
st2 key set notification_email YOUREMAIL
```
- `notification_email`: Email address that will receive notifcations.

## Actions

ref                                              | description
-------------------------------------------------|------------------------------------------------------------------------
Annotationtracks.annotationtracks                | Creates annotation trackfiles from the manual cluster data.
Annotationtracks.manual_clustering               | Cluster the given structural variants using manual coding based on similarity and functions.
Annotationtracks.hdbscan_clustering              | Cluster the given structural variants using HDBSCAN based on similarity and functions.
Annotationtracks.gene_data                       | Retrieve and filter gene data based on MANE status.
Annotationtracks.regulatory_data                 | Retrieve and filter regulatory data.
Annotationtracks.public_variants                 | Retrieve and filter known pathogenic/common structural variants.
Annotationtracks.dosage_sensitivity              | Retrieve and filter known dosage sensitive genes.
Annotationtracks.disease_related_variants        | Retrieve and filter litterature based disease related variants.
Annotationtracks.clinical_significance           | Retrieve and filter clinical significant structural variants.
Annotationtracks.acmg_genes                      | Retrieve and filter ACMG genes.                      
Annotationtracks.upload_tracks_to_gens           | Uploads annotation tracks to Gens.



## Workflows

ref                                                   | description
------------------------------------------------------|------------------------------------------------------------------------
Annotationtracks.update_structuralvariant_tracks_gens | Workflow for generating and updating structural variant annotationtracks in Gens.
Annotationtracks.update_genomic_elements              | Workflow for generating and updating genomic element annotationtracks in Gens.
Annotationtracks.update_informative_tracks_gens       | Workflow for generating and updating informative annotationtracks in Gens.

## Rules

ref                                                    | description
-------------------------------------------------------|------------------------------------------------------------------------
Annotationtracks.update_structuralvariants_tracks      | Generate a timer to update the annotationtracks for structural variants every friday at 18.
Annotationtracks.update_informative_tracks             | Generate a timer to update the informative (disease related SVs, clinical significant SVs, etc.) once every month


