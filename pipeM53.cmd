#!/bin/bash

# Pipeline from photometry to light curve statistics for M53J
# Ete 2017 - 1396 LAL, F. Habibi

# Note: The images should be previously regularised for mask

echo ''
echo ''
echo 'Photometric pipe line for OSER-M53'
echo ''


#RegDir=M53JReg_1
RegDir=M53JReg_2

CatDir=Catalogues

rm -r $RegDir/$CatDir
mkdir $RegDir/$CatDir/
cp *.cc *.h $RegDir/$CatDir/.

# Number of images for a given target or night
#imagenum=37 # first night
imagenum=20 # second night
echo "Number of images: " $imagenum
echo "  "

g++ makelist.cc -o ml
./ml $imagenum
cp imnum.list $RegDir/$CatDir/.


echo " Aperture photometry by SExtractor ... "
nhead=7 # number of parameters (and lines of meta data) to be in the catalgues
echo 'Number of parameters in default.param:' $nhead


#Making the reference catalog from the first image of the first night
assoclist=refM53J.cat #assocList.cat
sex M53JReg_1/defmaskReg0.fits -CATALOG_NAME $assoclist

g++ starCount.cc -o stcount
./stcount $nhead $assoclist
Nassoc=$( cat nassoc.dat )

###########################################
#mean back ground emission
mbg=9000

# margins for the images
marge=30

make clean
make
# making a reference synthetic image
rm synRef.fits
./synRef $Nassoc $assoclist $nhead $mbg $marge
sex synRef.fits,$RegDir/defmaskReg0.fits -CATALOG_NAME $assoclist
# reference image for the second night
#rm synRef.fits
#./synRef2 $Nassoc $assoclist $nhead $mbg $decalX $decalY
###########################################

cp $assoclist $RegDir/$CatDir/.


###########################################
# counting the number of detected stars (from the synthetiv reference image) saved in $assoclist
g++ starCount.cc -o stcount
./stcount $nhead $assoclist
Nassoc=$( cat nassoc.dat )
###########################################


#decalX=0 #shift in (pixel) x direction between the reference image and first image of the first nigh
decalX=-20.5 # same as above but for the second night
#decalY=0 #shift in (pixel) y direction between the reference image and first image of the first nigh
decalY=23 # same as above but for the second night

stepX=-0.054 # pixel shift in x direction for consecutive images (same fo both nights)
stepY=0.03 # first night
#stepY=0 # second nihgt


for i in $( cat imnum.list )
do
    rm  decalSynRef.fits
    # aligning the reference image with science image
    #./decalSynRef 0 0 -$stepX -$stepY $i
    #sex decalSynRef.fits,$RegDir/defmaskReg$i.fits -CATALOG_NAME $RegDir/$CatDir/im$i.cat
    ./decalSynRef decalX decalY $stepX $stepY $i $RegDir/defmaskReg$i.fits
    # detecting on the reference image and photometry on science images
    #sex M53JReg_1/defmaskReg0.fits,decalSynRef.fits -CATALOG_NAME $RegDir/$CatDir/im$i.cat
    sex synRef.fits,decalSienceIm.fits -CATALOG_NAME $RegDir/$CatDir/im$i.cat
done

echo  ''
echo 'Changing to directory' $RegDir/$CatDir/
cd $RegDir/$CatDir/


#echo 'Combining the catalogues ...'
rm manage
g++ manageCat1.cc -o manage
for i in $( cat imnum.list )
do
    ./manage im$i.cat imFidAssoc$i.cat $i $Nassoc $nhead $stepX $stepY $assoclist
done

echo 'Flux calibration and writting all measurements in allmeas.cat ... '
rm append
rm allmeas.cat
g++ append.cc -o append
for i in $( cat imnum.list )
do
    ./append imFidAssoc$i.cat $i $nhead $Nassoc
done

echo 'Computing the precison curve ... '
g++ computeLcStat.cc stat.cc -o cls
./cls $imagenum $Nassoc allmeas.cat LcStat.info


# example of spectral density computation for a star
#g++ computeAutoCorr.cc stat.cc -o cac -I ~/fink-0.39.3/include -L ~/fink-0.39.3/lib -lfftw3 -lm
#./cac corr343.dat 343 $Nassoc $imagenum

echo 'Computing the spectral density and autocorrelation for all light curves ...'
g++ computeAllAutoCorr.cc stat.cc -o caac -I ~/fink-0.39.3/include -L ~/fink-0.39.3/lib -lfftw3 -lm
./caac $imagenum $Nassoc allmeas.cat allCorr.dat

echo 'Computing the mean autocorr and spectrum for stable and variable stas
g++ computeMeanCorr.cc stat.cc -o cmc -I ~/fink-0.39.3/include -L ~/fink-0.39.3/lib -lfftw3 -lm
./cmc $imagenum $Nassoc


