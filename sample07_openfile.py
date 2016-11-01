#!/usr/bin/env python
# -*- coding: utf-8 -*-

import platform
if platform.python_version_tuple()[0] != "2":
    print "Sorry, this code runs on python 2"

import os, sys, codecs

# ファイルの連続処理の例
# import glob
# for filename in glob.glob( dataFolder + "/*.*"):
#	print filename
#	process_file( filename )

# ファイルを開く例（特にutf8の場合など）

for line in codecs.open('file_in.csv', mode='r', encoding='utf-8'):
	print line

# こういう書き方もあります
# with codecs.open('file_in.csv', mode='r', encoding='utf-8') as file_in:
#   for line in file_in:


# 書き出す例
# file_out = codecs.open('out.csv', 'w', 'utf-8')
# file_out.write( u"time,lat,lon,depth,mag\n" )

