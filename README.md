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

## Workflows

ref                                                  | description
-----------------------------------------------------|------------------------------------------------------------------------
st2_annotationtracks.update_annotationtracks_in_gens | Workflow for generating and updating annotationtrack files in Gens.

## Rules

ref                                                  | description
-----------------------------------------------------|------------------------------------------------------------------------
st2_annotationtracks.update_annotation_tracks        | Generate a timer to update the annotationtracks every friday at 18.


