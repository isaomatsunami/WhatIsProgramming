#!/usr/bin/env python
# -*- coding: utf-8 -*-

def isPrime(n):
  # 2以上n未満で割り切れたらFalseを返す
  for i in range(2, n):
    if n % i == 0:
      return False
  return True

a = 14414443
if isPrime(a):
  print( str(a) + "は素数です" )
else:
  print( str(a) + "は素数ではありません" )
