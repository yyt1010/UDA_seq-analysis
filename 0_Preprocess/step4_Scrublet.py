import scanpy as sc
import scrublet as scr
import numpy as np
anno=sc.read_h5ad("seurat2h5.h5ad") ##seurat2h5.h5ad represented a converted h5 file from seurat object merged_gene_mat.
scrub = scr.Scrublet(anno.X)
doublet_scores, predicted_doublets = scrub.scrub_doublets()
print(predicted_doublets)
print(len(predicted_doublets))
barcode=np.array(anno.obs_names)
doublet=np.array(predicted_doublets)
np.savetxt("scru_barcode.txt",barcode,fmt="%s")
np.savetxt("scru_doublet.txt",doublet,fmt="%s")

