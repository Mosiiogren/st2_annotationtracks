# StackStorm pack for creating and uploading annotationtrackfiles to Gens

## Installation
```bash
st2 pack install https://github.com/Mosiiogren/st2_annotationtracks
st2 pack config st2_annotationtracks
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

ref                                                  | description
-----------------------------------------------------|------------------------------------------------------------------------
st2_annotationtracks.annotationtracks                | Creates annotation trackfiles from the cluster data.
st2_annotationtracks.clustering                      | Cluster the given structural variants based on similarity and functions.
st2_annotationtracks.gene_data                       | Retrieve and filter gene data based on MANE status.
st2_annotationtracks.regulatory_data                 | Retrieve and filter regulatory data.
st2_annotationtracks.public_variants                 | Retrieve and filter known pathogenic/common structural variants.
st2_annotationtracks.dosage_sensitivity              | Retrieve and filter known dosage sensitive genes.
st2_annotationtracks.disease_related_variants        | Retrieve and filter litterature based disease related variants.
st2_annotationtracks.clinical_significance           | Retrieve and filter clinical isgnificant structural variants

## Workflows

ref                                                       | description
----------------------------------------------------------|------------------------------------------------------------------------
st2_annotationtracks.update_structuralvariant_tracks_gens | Workflow for generating and updating structural variant annotationtracks in Gens.
st2_annotationtracks.update_genomic_elements              | Workflow for generating and updating genomic element annotationtracks in Gens.
st2_annotationtracks.update_informative_tracks_gens       | Workflow for generating and updating informative annotationtracks in Gens.

## Rules

ref                                                        | description
-----------------------------------------------------------|------------------------------------------------------------------------
st2_annotationtracks.update_annotation_tracks              | Generate a timer to update the annotationtracks every friday at 18.


