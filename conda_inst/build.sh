#!/bin/bash
BINDIR="$PREFIX/bin"
#mkdir $BINDIR
cp -R $SRC_DIR/shamanapp $BINDIR
cp -R $SRC_DIR/KronaRShy $BINDIR
cp -R $SRC_DIR/shaman_bioblend $BINDIR
cp $BINDIR/shamanapp/www/shaman_conda.R $BINDIR/
chmod +x $BINDIR/shaman_conda.R
