mamba create -y -n genespace mamba openjdk orthofinder \
r-essentials r-devtools \
bioconductor-biostrings=2.66.0 bioconductor-rtracklayer=1.58.0

mkdir -p ./tools
cd ./tools
git clone https://github.com/wyp1125/MCScanX.git
cd MCScanX
make
export PATH=$PATH:/home/jbrand/tools/MCScanX
cd ..
conda activate genespace

# we were using this version 1.0.8 which can be installed like this
devtools::install_github("jtlovell/GENESPACE", upgrade = F, ref = 'd735cc9')
