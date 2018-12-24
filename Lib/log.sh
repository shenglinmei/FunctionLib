ls |while read aa; do if [ -d $aa ]; then echo $aa; fi; done


cat file2.txt |while read aa; do arr=(${aa//// });b=${arr[6]};echo $b ; done

 rsync -av -e ssh --exclude='*.rds' /home/meisl/Workplace/BMME/a.data/old/cellType/ meisl@pklab.med.harvard.edu://home/meisl/public_html/BMME/report/cellType/

ls ~/Workplace/BMME/a.data/conos/ada/conos/all_noBP/|grep .cell.cluster.ratio.png|while read aa; do arr=(${aa/.cell.cluster/ });b=${arr[0]}; echo bash /home/meisl/Workplace/BMME/a.data/conos/ada/html.sh $b /home/meisl/Workplace/BMME/a.data/conos/ada/conos/all_noBP/ /home/meisl/Workplace/BMME/a.data/conos/ada/conos_tumor/html/ALL_samples/$b".html" >> sh.sh ;done



PATH=/home/meisl/.local/bin/:$PATH
Rscript -e "rmarkdown::render('r.rmd',output_file='output.html')"

rsync -av -e ssh --exclude='*.rds' --exclude '*bin' ./Tumor_all_new meisl@mendel.med.harvard.edu:/home/meisl/Workplace/BMME/a.data/conos/ada/

Rscript -e  "rmarkdown::render('r.rmd',params = list(appname   = 'Whole_Tumor_new_t2_PCA_15_30_500_conos.rds',path = '/home/meisl/Workplace/BMME/a.data/conos/ada/conos_tumor/Whole_Tumor_new/'),output_file='output.html')"

