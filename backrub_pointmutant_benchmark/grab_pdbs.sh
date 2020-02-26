rsync -av --include=*/*/$1/*** \
    --include=*/ \
    --exclude=* \
    --include=pdbs_*.tar.gz \
    ckrivacic@guybrush.ucsf.edu:/kortemmelab/home/ckrivacic/intelligent_design/sidechain_directed_design_pyrosetta/backrub_pointmutant_benchmark/output/ \
    output/
