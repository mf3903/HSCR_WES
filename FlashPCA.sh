project=""

SUFFIX=WES

# Run FlashPCA2 on the pruned genotypes
dx run swiss-army-knife \
    -iin="${project}:/Genotypes/Rare Variant PCs/ukb23159_exome_filtered_for_fine_structure_PCs.${SUFFIX}.bed" \
    -iin="${project}:/Genotypes/Rare Variant PCs/ukb23159_exome_filtered_for_fine_structure_PCs.${SUFFIX}.bim" \
    -iin="${project}:/Genotypes/Rare Variant PCs/ukb23159_exome_filtered_for_fine_structure_PCs.${SUFFIX}.fam" \
    -iimage_file="${project}:/Software/flashpca2.tar.gz" \
    -imount_inputs=true \
    -icmd="
    flashpca \
        --bfile ukb23159_exome_filtered_for_fine_structure_PCs.${SUFFIX} \
        --numthreads 16 \
        --ndim 20 \
        --outpc ukb23159_exome_filtered_for_fine_structure_PCs.pcs_flashpca.${SUFFIX}.txt \
        --outvec ukb23159_exome_filtered_for_fine_structure_PCs.eigenvectors_flashpca.${SUFFIX}.txt \
        --outval ukb23159_exome_filtered_for_fine_structure_PCs.eigenvalues_flashpca.${SUFFIX}.txt \
        --outpve ukb23159_exome_filtered_for_fine_structure_PCs.pve_flashpca.${SUFFIX}.txt
    " \
    --destination "${project}:/Genotypes/Rare Variant PCs/" \
    --instance-type mem1_ssd1_v2_x16 \
    --priority low \
    --name "Estimate Fine Structure PCs" \
    --yes
