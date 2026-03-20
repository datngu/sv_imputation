
## chain file
mkdir -p data/genomes


## chain file
wget https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/GRCh37_to_GRCh38.chain.gz -O data/genomes/GRCh37_to_GRCh38.chain.gz

## hg19 and hg38 reference genomes
wget -qO- http://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz | gunzip -c > data/genomes/GRCh37.fa


wget -qO- http://ftp.ensembl.org/pub/release-106/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | gunzip -c > data/genomes/GRCh38.fa


## index genomes
samtools faidx data/genomes/GRCh37.fa
samtools faidx data/genomes/GRCh38.fa

