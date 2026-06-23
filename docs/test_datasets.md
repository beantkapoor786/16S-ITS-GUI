# Test Datasets for RAPID-16S

Two datasets are provided below. A quick test dataset for first-time users, and a real soil dataset for users who want to test with field microbiome data.

---

## Option 1 - Quick Test Dataset (Recommended for First-Time Users)

**DADA2 MiSeq SOP Dataset** — Mouse gut, 16S V4 region, 20 samples, 2×250 bp Illumina MiSeq  
This is the official test dataset from the DADA2 authors, widely used in tutorials and benchmarking. It is small (~40 MB), well-characterized, and runs through the full RAPID-16S pipeline in minutes.

- **Primers:** 515F / 806R (V4 region)
- **Read length:** 2×250 bp paired-end
- **Samples:** 20 mouse gut samples
- **Paper:** Kozich et al. (2013) *Appl. Environ. Microbiol.* [doi:10.1128/AEM.01043-13](https://doi.org/10.1128/AEM.01043-13)

### Download

**wget:**
```bash
wget http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
unzip MiSeqSOPData.zip
```

**curl:**
```bash
curl -O http://www.mothur.org/w/images/d/d6/MiSeqSOPData.zip
unzip MiSeqSOPData.zip
```

Point RAPID-16S to the `MiSeq_SOP/` folder after unzipping.

### Suggested RAPID settings for this dataset
| Parameter | Value |
|---|---|
| Forward pattern | `_R1_` |
| Reverse pattern | `_R2_` |
| truncLen (F) | 240 |
| truncLen (R) | 160 |
| Taxonomy database | SILVA v138.1 |

---

## Option 2 — Real Soil Dataset

**Differently Managed Soils — Comparative primer study** — Soil, 16S V4 region, 515F-806R, Illumina MiSeq 2×150 bp  
Agricultural soils collected from differently managed fields. Data was originally generated to compare primer pairs for prokaryotic and archaeal nitrifier detection in soils, processed with DADA2 via QIIME 2.

- **Primers:** 515F / 806R (V4 region)
- **Read length:** 2×150 bp paired-end
- **Sample type:** Agricultural soil
- **NCBI BioProject:** [PRJNA831877](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA831877)
- **ENA mirror:** [https://www.ebi.ac.uk/ena/browser/view/PRJNA831877](https://www.ebi.ac.uk/ena/browser/view/PRJNA831877)
- **Paper:** Duan et al. (2023) *Front. Microbiol.* [doi:10.3389/fmicb.2023.1140487](https://doi.org/10.3389/fmicb.2023.1140487)

### Download Options

**Option A — ENA browser (easiest, no terminal needed)**

Go to [https://www.ebi.ac.uk/ena/browser/view/PRJNA831877](https://www.ebi.ac.uk/ena/browser/view/PRJNA831877), click **Show Column Selection**, enable **FASTQ files (FTP)**, then click the download links.

---

**Option B — wget**

```bash
# Step 1: Download the file report from ENA
wget -O filereport.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA831877&result=read_run&fields=run_accession,fastq_ftp&format=tsv&download=true"

# Step 2: Extract FASTQ URLs and download all files
cut -f2 filereport.tsv | tail -n +2 | tr ';' '\n' | sed 's|^|ftp://|' | wget -i -
```

---

**Option C — curl**

```bash
# Step 1: Download the file report
curl -o filereport.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA831877&result=read_run&fields=run_accession,fastq_ftp&format=tsv&download=true"

# Step 2: Download each FASTQ file
cut -f2 filereport.tsv | tail -n +2 | tr ';' '\n' | while read url; do
    curl -O "ftp://$url"
done
```

---

**Option D — SRA Toolkit (NCBI)**

```bash
# Install SRA Toolkit: https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit

# Get all accession IDs for this BioProject
esearch -db sra -query PRJNA831877 | efetch -format runinfo | \
  cut -d',' -f1 | tail -n +2 > accessions.txt

# Download and convert to paired-end FASTQ
prefetch --option-file accessions.txt
fasterq-dump --split-files SRR*
```

---

**Option E — ffq (Python, recommended for scripted workflows)**

```bash
# Install: pip install ffq
ffq --ftp PRJNA831877 | jq -r '.[] | .url' | xargs -n 1 wget
```

---

**Option F — Aspera (fastest for large transfers)**

```bash
# Install IBM Aspera Connect: https://www.ibm.com/aspera/connect/
# ENA Aspera URLs use era-fasp@ prefix — get them from the file report:
wget -O filereport_aspera.tsv "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA831877&result=read_run&fields=run_accession,fastq_aspera&format=tsv&download=true"

cut -f2 filereport_aspera.tsv | tail -n +2 | tr ';' '\n' | while read url; do
    ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
      era-fasp@$url .
done
```

### Suggested RAPID settings for this dataset
| Parameter | Value |
|---|---|
| Forward pattern | `_1.fastq` |
| Reverse pattern | `_2.fastq` |
| truncLen (F) | 150 |
| truncLen (R) | 140 |
| Taxonomy database | SILVA v138.1 |

> **Note:** Because this dataset uses 2×150 bp reads (shorter than 2×250 bp), set `truncLen` conservatively — inspect the quality profiles in Step 2 before proceeding.
