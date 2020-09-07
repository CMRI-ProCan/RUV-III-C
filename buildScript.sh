R CMD build .
R CMD check RUVIIIC_1.0.15.tar.gz --as-cran
mkdir repack/
cd repack/
tar -zxf ../RUVIIIC_1.0.16.tar.gz
rm -rf ./RUVIIIC/.circleci ./RUVIIIC/version ./RUVIIIC/docker ./RUVIIIC/buildScript.sh ./RUVIIIC/README.md ./RUVIIIC/repack
tar -czf ./RUVIIIC_1.0.16.tar.gz RUVIIIC
