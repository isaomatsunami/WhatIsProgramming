#!/usr/bin/env python
# -*- coding: utf-8 -*-

def printSqrt(n):
  print( n * n )
  return n*n
job = printSqrt
items = [3, 6, 2, -1]
res = map( printSqrt, items )
print( list(res) )
